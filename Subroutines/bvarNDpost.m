% Written by:
% Uriel Braham
% uriel.braham@googlemail.com

function [result, VAR] = bvarNDpost(VARoption,BVARoption,result,VAR,DATA)

T = VAR.T;
M = VAR.M;
N = VAR.N;
K = VAR.K;
k = VAR.k;
q = VAR.q;
p = VARoption.p;
VAR.k = M*K;


Y = VAR.Yreg;
X = VAR.Xreg;
Yexo = DATA.Yexo;
Yendo = DATA.Yendo;

A_draws     = result.A_draws;
SIGMA_draws = result.SIGMA_draws;

A_post             = squeeze(quantile(A_draws,0.5));
a_post             = vec(A_post);
A_SD               = squeeze(std(A_draws,1)); 
SIGMA_post         = squeeze(quantile(SIGMA_draws,0.5));

% SSE
SSE = (Y - X*A_post)'*(Y - X*A_post);

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

%==========================================================================
% RESULTS
%==========================================================================
result.A_post      = A_post;
result.a_post      = a_post;
result.A_SD        = A_SD;
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
VAR.K           = K;
VAR.k           = k;
VAR.T           = T;
VAR.Yreg        = Y;
VAR.Xreg        = X;

