% Written by:
% Uriel Braham
% uriel.braham@googlemail.com

function [result, VAR] = bvarDiffusepost(VARoption,DATA)

%**************************************************************************
%   FUNCTION TO ESTIMATE A LINEAR WITH DIFFUSE PRIOR
%--------------------------------------------------------------------------
%   INPUT: 
%   - Structure 'VARoption'
%   - Structure 'DATA'
%--------------------------------------------------------------------------
%   OUTPUT: 
%   - Structure 'result' with estimation result of the BAVR
%   - Structure 'VAR' with estimation information
%--------------------------------------------------------------------------
%
%**************************************************************************

%----------------------------CHECK INPUTS----------------------------------
if ~exist('VARoption')
    error('You need to provide VAR options in "VARoption".');
end
% Check they are not empty
Yendo       = DATA.Yendo;
Yexo        = DATA.Yexo;
if isempty(Yendo)
    error('You need to provide endogenous variables in "DATA".');
end
%--------------------------------------------------------------------------

%*************************************************************************%
%                UNPACK structures 'VARoption' and 'VAR'                  %
%*************************************************************************%
% OPTIONS
constant    = VARoption.constant;
p           = VARoption.p;
trend       = VARoption.trend;
nsave       = VARoption.nsave;
nburn       = VARoption.nburn;
ntot        = nburn + nsave;
it_print = 2000;

% Notes on dimensions:
% M=number of endogenous variables
% N=number of exogenous (deterministic) variables (without constant and/or trend) 
% K=q+M*p total number of coefficients in each equation i=1,...,M
% k=M*K total number of coefficients

% Get initial dimensions
[Traw, M] = size(Yendo);
N = size(Yexo,2);
q = constant + trend + N;


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
k = M*K;
% 
% Matricies for estimation
Y = Y1;
X = X1;

%*************************************************************************%
%                             POSTERIORS                                  %
%*************************************************************************%
% Posterior mean of coefficients:
A_OLS = (X'*X)\(X'*Y);
A_post = A_OLS;
a_OLS = vec(A_OLS);
a_post = a_OLS;
SSE = (Y - X*A_post)'*(Y - X*A_post);
SIGMA_OLS = SSE./(T-K);
SIGMA_post = SIGMA_OLS;
a_COV = (kron(SIGMA_post, inv(X'*X)));
a_SD  =diag(sqrt(a_COV)); 

% R-SQUARED
muY = mean(Y)';
SSE_d = diag((Y-X*A_post)'*(Y-X*A_post));
R2 = 1 - diag(SSE_d)./diag((Y - muY')'*(Y - muY'));
R2 = diag(R2);

% ROOT MEAN SQUARED ERROR
RMSE = sqrt(diag(SIGMA_post));

% FIT
SS = 1/(T-p+1)*(Y-X*A_post)'*(Y-X*A_post);
FIT = 1 - diag(SS)./diag(cov(diff(Yendo)));

% FITTED VALUES
Yendo_fitted = X*A_post;

% STABILITY
[S, maxeig] = stability(A_post,M,p,q);

%*************************************************************************%
%                             RESULTS                                     %
%*************************************************************************%
% STRUCTURE OF RESULTS
result.A_post      = A_post;
result.a_post      = a_post;
result.a_SD        = a_SD;
result.SIGMA_post  = SIGMA_post;
result.SSE         = SSE;
result.RMSE        = RMSE;
result.R2          = R2;
result.FIT         = FIT;
%result.IC          = IC;
result.FITTED      = Yendo_fitted;
result.S           = S;
result.maxeig      = maxeig;

% STRUCTURE OF DATA
VAR.M           = M;
VAR.N           = N;
VAR.q           = q;
VAR.K           = K;
VAR.k           = k;
VAR.T           = T;
VAR.Yreg        = Y;
VAR.Xreg        = X;
%--------------------------------------------------------------------------
