% Written by:
% Uriel Braham
% Univerity of Bayreuth
% Universitätsstraße 30
% 95445 Bayreuth
% uriel.braham@googlemail.com

function [result, VAR] = varOLS(VARoption,DATA);

%**************************************************************************
%   FUNCTION TO ESTIMATE A LINEAR WITH OLS
%--------------------------------------------------------------------------
%   INPUT: Structure 'DATA' with
%   - Yendo: (T x M) matrix of endogenous variables
%   - Yexo: (T x N) matrix of exogenous variables
%          Structure 'VARoption' with
%   - constant: 1=constant, 0=no constant
%   - trend:  1=linear trend, 0=no linear trend
%   - p: lag-order
%--------------------------------------------------------------------------
%   OUTPUT: STRUCTURE 'OLS'
%   - A_OLS: (K x M) matrix of coefficients
%   - A_SE: (K x M) matrix of standard errors
%   - SIGMA_OLS: (M x M) matrix of estimated resiudal variances
%   - R2: (M x 1) vector of R-Squares
%   - FIT: (M x 1) vector of relative in-sample fit
%   - IC: (3 x 1) vector of Akaike, Schwarz and Hannan-Quinn Info. Criteria
%--------------------------------------------------------------------------
%
%**************************************************************************

%----------------------------DATA PREPERATION------------------------------
% OPTIONS
constant    = VARoption.constant;
p           = VARoption.p;
trend       = VARoption.trend;
% DATA
Yendo       = DATA.Yendo;
Yexo        = DATA.Yexo;

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
% 
% Matricies for estimation
Y = Y1;
X = X1;

%-------------------------------ESTIMATION---------------------------------

% OLS ESTIMATES
A_OLS = (X'*X)\(X'*Y);              
SSE = (Y - X*A_OLS)'*(Y - X*A_OLS);
SIGMA_OLS = SSE./(T-p*M-q);
% ML ESTMATES
%A_MLE = (X'*X)\(X'*Y);              
%SSE = (Y - X*A_MLE)'*(Y - X*A_MLE);
%SIGMA_MLE = SSE./T;

% OLS STANDARD ERRORS
a_SE = kron(SIGMA_OLS, inv((X'*X)));
A_SE = reshape(diag(a_SE),K,M); 
A_SE = sqrt(A_SE);
%a_SE = kron(SIGMA_MLE, (X'*X)\eye(size(X,1)));
%A_SE = reshape(diag(a_SE),K,M); 

% R-SQUARED
muY = mean(Y)';
SSE_d = diag((Y-X*A_OLS)'*(Y-X*A_OLS));
R2 = 1 - diag(SSE_d)./diag((Y - muY')'*(Y - muY'));
R2 = diag(R2);

% ROOT MEAN SQUARED ERROR
RMSE = sqrt(diag(SIGMA_OLS));

% FIT
SS = 1/(T-p+1)*(Y-X*A_OLS)'*(Y-X*A_OLS);
FIT = 1 - diag(SS)./diag(cov(diff(Yendo)));

% FITTED VALUES
Yendo_fitted = X*A_OLS;

% STABILITY
[S, maxeig] = stability(A_post,M,p,q);

% INFORMATION CRITERIA: 
% Note: IC have to be evaluated at the MLE (see Kilian and Lütkepohl 2017,
% p.56).
k = size(vec(A_OLS),1);
SIGMA_MLE = SSE./T; 
AICcrit = log(det(SIGMA_MLE))+k*2/(T);    	   	    % AIC value
HQCcrit = log(det(SIGMA_MLE))+k*2*log(log(T))/(T);  % HQC value
SICcrit = log(det(SIGMA_MLE))+k*log(T)/(T);         % SIC value
IC = [AICcrit; HQCcrit; SICcrit];

%--------------------------------RESULTS-----------------------------------
% STRUCTURE OF OLS RESULTS
result.A_OLS       = A_OLS;
result.A_SE        = A_SE;
result.SIGMA_OLS   = SIGMA_OLS;
result.RMSE    = RMSE;
result.R2      = R2;
result.FIT     = FIT;
result.IC      = IC;
result.FITTED  = Yendo_fitted;
result.S           = S;
result.maxeig      = maxeig;

% STRUCTURE OF DATA
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


