% Written by:
% Uriel Braham
% uriel.braham@googlemail.com

function [result, VAR] = bvarNDGibbs(VARoption,BVARoption,VAR,DATA)

%**************************************************************************
%   FUNCTION TO DRAW FROM POSTERIOR DISTIBUTIONS WITH GIBBS SAMPLER
%--------------------------------------------------------------------------
%   INPUT: 
%   - Structure 'VARoption'
%   - Structure 'BVARoption'
%   - Structure 'VAR'
%   - Structure 'DATA'
%--------------------------------------------------------------------------
%   OUTPUT: 
%   - Structure 'result' with estimation results
%   - Structure 'VAR' with VAR information 
%--------------------------------------------------------------------------
%
%**************************************************************************

%*************************************************************************%
%                            CHECK INPUTS                                 %
%*************************************************************************%
if ~exist('VARoption')
    error('You need to provide VAR options in "VARoption".');
end
if ~exist('BVARoption')
    error('You need to provide BVAR options in "BVARoption".');
end

% Check they are not empty
Yendo       = DATA.Yendo;
Yexo        = DATA.Yexo;
if isempty(Yendo)
    error('You need to provide endogenous and exogenous variables variables in "DATA".');
end

%*************************************************************************%
%                UNPACK structures 'VARoption' and 'VAR'                  %
%*************************************************************************%
constant        = VARoption.constant;
p               = VARoption.p;
trend           = VARoption.trend;
nsave           = VARoption.nsave;
nburn           = VARoption.nburn;
ntot            = nburn + nsave;

%*************************************************************************%
%                           DATA PREPERATIONS                             %
%*************************************************************************%
% Notes on dimensions:
% M=number of endogenous variables
% N=number of exogenous (deterministic) variables (without constant and/or trend) 
% K=q+M*p total number of coefficients in each equation i=1,...,M
% k=M*K total number of coefficients

% Get initial dimensions
[Traw, M] = size(Yendo);
N = size(Yexo,2);
q = constant + trend + N;
K = M*p + q;
k = M*K;

Y1 = Yendo;
Y2 = Yendo;

% Generate lagged Y matrix. This will be part of the X matrix
Ylag=zeros(Traw,M*p);
for ii=1:p
    Ylag(p+1:Traw,(M*(ii-1)+1):M*ii)=Y2(p+1-ii:Traw-ii,1:M);
end

% Generate the regressor matrix
if constant ==1
    if isempty(Yexo)==1
        if trend==1
                X1 = [Ylag(p+1:Traw,:) ones(Traw-p,1) (p+1:Traw)'];
        elseif trend==0 
                X1 = [Ylag(p+1:Traw,:) ones(Traw-p,1) ];  
         end 
    elseif isempty(Yexo)==0
         if trend ==1
                X1 = [Ylag(p+1:Traw,:) ones(Traw-p,1) ((p+1):Traw)' Yexo(p+1:Traw,:)];
         elseif trend==0 
                X1 = [Ylag(p+1:Traw,:) ones(Traw-p,1) Yexo(p+1:Traw,:) ]; 
         end
    end 
elseif constant==0    
    if isempty(Yexo)==1
        if trend==1
                X1 = [Ylag(p+1:Traw,:) (p+1:Traw)' ];
        elseif trend==0 
                X1 = [Ylag(p+1:Traw,:)];  
        end 
    elseif isempty(Yexo)==0
        if trend ==1
                X1 = [Ylag(p+1:Traw,:) ((p+1):Traw)' Yexo(p+1:Traw,:)];
        elseif trend==0 
                X1 = [Ylag(p+1:Traw,:) Yexo(p+1:Traw,:)]; 
        end
    end  
end


% Delete first p-rows to match the dimensions of X matrix
Y1 = Y1(p+1:Traw,:);
% Traw number of obs. in initial data. T is the number of actual 
% time series observations of Y and X. 
T = Traw - p;
% Get size of final matrix X
[Traw3 K] = size(X1); 
% 
% Matricies for estimation
Y = Y1;
X = X1;


%*************************************************************************%
%                           OLS QUANTITIES                                %
%*************************************************************************%
A_OLS = (X'*X)\(X'*Y);
a_OLS = vec(A_OLS);
SSE_OLS = (Y - X*A_OLS)'*(Y - X*A_OLS);   
SIGMA_OLS = SSE_OLS./(T-K);

%*************************************************************************%
%                              PRIOR                                      %
%*************************************************************************%
[NDprior] = bvarNDprior(VARoption,BVARoption,DATA);
A_prior     = NDprior.A_prior;
a_prior     = NDprior.a_prior;
v_prior     = NDprior.v_prior;
PHI_prior   = NDprior.PHI_prior;
Yd          = NDprior.Yd;
Xd          = NDprior.Xd;
Td          = NDprior.Td;

%*************************************************************************%
%                             POSTERIORS                                  %
%*************************************************************************%
% Append dummy observations
Ys           = [Y; Yd];
Xs           = [X; Xd];
Ts           = T + Td;
invPHI_prior=diag(1./diag(PHI_prior));
v_post       = T + Td + v_prior;
SCALE_post=(Ys-Xs*A_OLS)'*(Ys-Xs*A_OLS);
% Correct potential asymmetries due to rounding errors from Matlab
C=chol(nspd(SCALE_post));
SCALE_post=C'*C;

%*************************************************************************%
%                             GIBBS SAMPLER                               %
%*************************************************************************%
% Storage space for posterior draws
a_draws         = zeros(nsave,K*M);   
A_draws         = zeros(nsave,K,M);  
SIGMA_draws     = zeros(nsave,M,M);

% Initial Draw
a_draw          = a_OLS;
A_draw          = A_OLS;
disp('')
disp('Posterior simulation for BVAR coefficients (Gibbs sampler):')
tic;
%Start Gibbs "loop"
for irep = 1:ntot  % loop over irep
    
    if mod(irep,1000) == 0 % print iterations
        disp([num2str(irep) ' / ' num2str(ntot) ' draws']);
        %toc;
    end
    
    %======================== POSTERIOR DRAWS =============================
    % 1.) Draw SIGMA|a,y ~ iW(SCALE_post,v_post)
    SCALE_post=(Ys-Xs*A_draw)'*(Ys-Xs*A_draw);
    % Correct potential asymmetries due to rounding errors from Matlab
    C=chol(nspd(SCALE_post));
    SCALE_post=C'*C;
    SIGMA_draw=IW(SCALE_post,v_post);

    % 2.) Draw a|SIGMA,y ~ N(a_post, PHI_post)
    C=chol(nspd(SIGMA_draw));
    invC=C\speye(M);
    invSIGMA_draw=invC*invC';

    invPHI_post=invPHI_prior+kron(invSIGMA_draw,Xs'*Xs);
    C=chol(nspd(invPHI_post));
    invC=C\speye(k);
    PHI_post=invC*invC';

    a_post=PHI_post*(invPHI_prior*a_prior+kron(invSIGMA_draw,Xs')*vec(Ys));
    a_draw=a_post+chol(nspd(PHI_post),'lower')*mvnrnd(zeros(k,1),eye(k))';
    A_draw = reshape(a_draw,K,M);
    
    
    %======================================================================
    if irep > nburn        
    %Save draws of the parameters
    a_draws(irep-nburn,:)       = a_draw;
    A_draws(irep-nburn,:,:)     = A_draw;
    SIGMA_draws(irep-nburn,:,:) = SIGMA_draw; 
    
    % End saving  
    end
% End Gibbs loop
end

A_post             = squeeze(quantile(A_draws,0.5));
A_SD               = squeeze(std(A_draws,1)); 
SIGMA_post         = squeeze(quantile(SIGMA_draws,0.5));

%*************************************************************************%
%                             RESULTS                                     %
%*************************************************************************%
result.A_draws     = A_draws;
result.a_draws     = a_draws;
result.SIGMA_draws = SIGMA_draws;
result.A_postSD    = A_SD;
result.SIGMA_post  = SIGMA_post;
VAR.M           = M;
VAR.N           = N;
VAR.q           = q;
VAR.k           = k;
VAR.K           = K;
VAR.Ys          = Ys;
VAR.Xs          = Xs;
VAR.Ts          = Ts;
VAR.Yd          = Yd;
VAR.Xd          = Xd;
VAR.Td          = Td;
VAR.Yreg        = Y;
VAR.Xreg        = X;
VAR.T           = T;
