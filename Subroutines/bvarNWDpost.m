% Written by:
% Uriel Braham
% uriel.braham@googlemail.com


function [result, VAR] = bvarNWDpost(VARoption,BVARoption,VAR,DATA)

%**************************************************************************
%   FUNCTION TO ESTIMATE A LINEAR BVAR WITH MINNESOTA PRIOR AND DUMMY 
%   OBSERVATIONS
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
if ~exist('BVARoption')
    error('You need to provide BVAR options in "BVARoption".');
end
if ~exist('VAR')
    error('You need to provide VAR information in "VAR".');
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
il1         = 1/lambda1;
il2         = 1/lambda2;
il3         = 1/lambda3;
il4         = 1/lambda4;
il5         = 1/lambda5;
il6         = 1/lambda6;

%*************************************************************************%
%                           DATA PREPERATIONS                             %
%*************************************************************************%
% Notes on dimensions:
% M=number of endogenous variables
% N=number of exogenous (deterministic) variables (without constant and/or trend) 
% K=q+M*p total number of coefficients in each equation i=1,...,M
% k=M*K total number of coefficients

Yendo = DATA.Yendo;
Yexo = DATA.Yexo;
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
A_OLS       = (X'*X)\(X'*Y);        
SSE_OLS     = (Y - X*A_OLS)'*(Y - X*A_OLS);   
SIGMA_OLS   = SSE_OLS./(T-K);
muYendo     = mean(Yendo);

%*************************************************************************%
%                              DUMMIES                                    %
%*************************************************************************%
[NWDprior]  = bvarNWDprior(VARoption,BVARoption,DATA);
Yd          = NWDprior.Yd;
Xd          = NWDprior.Xd;
Td          = NWDprior.Td;

% Append the artifical data to the actual data
Ys          = [Y; Yd];
Xs          = [X; Xd];
Ts          = T + Td;

%*************************************************************************%
%                             POSTERIORS                                  %
%*************************************************************************%
% Conditional (posterior) mean of VAR coefficients
A_post      = (Xs'*Xs)\(Xs'*Ys);          
a_post      = vec(A_post);                
XXs         = Xs'*Xs;
invXXs      = XXs\eye(cols(XXs));        

% Posterior mean of Variance-Covariance matrix:
Us          = (Ys - Xs*A_post); 
SCALE_post  = (Us'*Us);
SIGMA_post  = (1/(T+Td- K - M - 1))*SCALE_post;
v_post      = Ts - K + 2; 

%*************************************************************************%
%                         FURTHER STATISTICS                              %
%*************************************************************************%

% SSE
SSE         = (Y - X*A_post)'*(Y - X*A_post);

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

% (LOG-) MARGINAL LIKELIHOOD
[lnmarglik] = bvarNWDmarglik(Xd,Yd,Xs,Ys);

% STABILITY
[S, maxeig] = stability(A_post,M,p,q);

%*************************************************************************%
%                             RESULTS                                     %
%*************************************************************************%
% STRUCTURE OF RESULTS
result.A_post      = A_post;
result.a_post      = a_post;
result.SCALE_post  = SCALE_post;
result.v_post      = v_post;
result.SIGMA_post  = SIGMA_post;
result.RMSE        = RMSE;
result.SSE         = SSE;
result.R2          = R2;
result.FIT         = FIT;
result.FITTED      = Yendo_fitted;
result.lnmarglik   = lnmarglik;
result.S           = S;
result.maxeig      = maxeig;
% STRUCTURE OF DATA
VAR.M           = M;
VAR.N           = N;
VAR.q           = q;
VAR.K           = K;
VAR.k           = k;
VAR.Yd          = Yd;
VAR.Xd          = Xd;
VAR.Td          = Td;
VAR.Ys          = Ys;
VAR.Xs          = Xs;
VAR.Ts          = Ts;
VAR.Yreg        = Y;
VAR.Xreg        = X;
VAR.T           = T;
%--------------------------------------------------------------------------
