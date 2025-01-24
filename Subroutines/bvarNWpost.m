% Written by:
% Uriel Braham
% uriel.braham@googlemail.com


function [result, VAR] = bvarNWpost(VARoption,BVARoption,VAR,DATA)

%----------------------------CHECK INPUTS----------------------------------
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
%--------------------------------------------------------------------------

%*************************************************************************%
%                UNPACK structures 'VARoption' and 'VAR'                  %
%*************************************************************************%
constant    = VARoption.constant;
p           = VARoption.p;
trend       = VARoption.trend;
lambda1     = BVARoption.lambda1;
lambda2     = BVARoption.lambda2;
lambda3     = BVARoption.lambda3;
lambda4     = BVARoption.lambda4;
lambda5     = BVARoption.lambda5; 
lambda6     = BVARoption.lambda6;
BVARoption.delta         = ones(length(BVARoption.stationary),1); 
BVARoption.delta(ismember(BVARoption.stationary, 'stationary')) = 0;  
deltai      = BVARoption.delta'; 

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
                X1 = [Ylag(p+1:Traw,:) ((p+1):Traw)' Yexo(p+1:Traw,:) ];
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
%                              OLS QUANTITIES                             %
%*************************************************************************%
A_OLS = (X'*X)\(X'*Y);        
SSE_OLS = (Y - X*A_OLS)'*(Y - X*A_OLS);   
SIGMA_OLS = SSE_OLS./(T-K);
muYendo = mean(Yendo);

%*************************************************************************%
%                               PRIORS                                    %
%*************************************************************************%
[NWprior] = bvarNWprior(VARoption,BVARoption,DATA);

A_prior     = NWprior.A_prior;
a_prior     = NWprior.a_prior;
v_prior     = NWprior.v_prior;
PHI_prior   = NWprior.PHI_prior;
SCALE_prior = NWprior.SCALE_prior;
SIGMA_prior = NWprior.SIGMA_prior;
Yd          = NWprior.Yd;
Xd          = NWprior.Xd;
Td          = NWprior.Td;

%*************************************************************************%
%                             POSTERIORS                                  %
%*************************************************************************%
% Append dummy observations
Ys              = [Y; Yd];
Xs              = [X; Xd];
Ts              = rows(Ys);

invPHI_prior    = diag(1./diag(PHI_prior));

invPHI_post     = invPHI_prior+Xs'*Xs;
C               = chol(nspd(invPHI_post),'Lower')';
invC            = C\speye(K);
PHI_post        = invC*invC';
%PHI_post = (invPHI_prior+X'*X)\speye(K);

A_post          = PHI_post*(invPHI_prior*A_prior+Xs'*Ys);
a_post          = vec(A_post);
    
v_post          = T + Td +v_prior;

SCALE_post      = Ys'*Ys+SCALE_prior+A_prior'*invPHI_prior*A_prior-A_post'*invPHI_post*A_post;
SCALE_post      = nspd(SCALE_post);
SIGMA_post      = SCALE_post./(v_post-M-1);

%*************************************************************************%
%                         FURTHER STATISTICS                              %
%*************************************************************************%

% SSE
SSE = (Y - X*A_post)'*(Y - X*A_post);

% Check for stability of the VAR model
[S]         =stability(A_post,M,p,q);

% ROOT MEAN SQUARED ERROR
RMSE        = sqrt(diag(SIGMA_post));

% R-SQUARED
muY         = mean(Y)';
SSE_d       = diag((Y-X*A_post)'*(Y-X*A_post));
R2          = 1 - diag(SSE_d)./diag((Y - muY')'*(Y - muY'));
R2          = diag(R2);

% FIT
SS          = 1/(T-p+1)*(Y-X*A_post)'*(Y-X*A_post);
FIT         = 1 - diag(SS)./diag(cov(diff(Yendo)));

% FITTED VALUES
Yendo_fitted= X*A_post;

% STABILITY
[S, maxeig] = stability(A_post,M,p,q);

%*************************************************************************%
%                             RESULTS                                     %
%*************************************************************************%
result.A_post      = A_post;
result.a_post      = a_post;
result.PHI_post    = PHI_post;
result.SCALE_post  = SCALE_post;
result.v_post      = v_post;
result.SIGMA_post  = SIGMA_post;
result.RMSE        = RMSE;
result.SSE         = SSE;
result.R2          = R2;
result.FIT         = FIT;
result.FITTED      = Yendo_fitted;
result.S           = S;
result.maxeig      = maxeig;

% STRUCTURE OF DATA
VAR.M           = M;
VAR.N           = N;
VAR.q           = q;
VAR.k           = k;
VAR.K           = K;
VAR.Yd          = Yd;
VAR.Xd          = Xd;
VAR.Td          = Td;
VAR.Ys          = Ys;
VAR.Xs          = Xs;
VAR.Ts          = Ts;
VAR.Yreg        = Y;
VAR.Xreg        = X;
VAR.Y           = Y;
VAR.X           = X;
VAR.T           = T;
%--------------------------------------------------------------------------

