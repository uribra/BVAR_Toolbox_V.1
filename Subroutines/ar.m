% Written by:
% Uriel Braham
% Univerity of Bayreuth
% Universit‰tsstraﬂe 30
% 95445 Bayreuth
% uriel.braham@googlemail.com

function [a_OLS]=ar(Yendo,p);

%**************************************************************************
%   FUNCTION TO ESTIMATE (UNIVARIATE) AR(p) MODEL
%--------------------------------------------------------------------------
%   INPUT: 
%   - Yendo: (T x 1) vector of endogenous variable
%   - p: lag-order
%--------------------------------------------------------------------------
%   OUTPUT: 
%   - A_OLS: (K x M) matrix of coefficients
%   - A_SE: (K x M) matrix of standard errors
%   - SIGMA_OLS: (M x M) matrix of estimated resiudal variances
%--------------------------------------------------------------------------
%
%**************************************************************************

sigma_sq    = zeros(M,1);            % vector to store residual variances
delta       = zeros(M,1);
%AR(p) model for each equation i=1,2,...,M 
for i = 1:M
    ytemp=Yendo(:,i);
    xtemp=[];
    for j=1:p
    xtemp=[xtemp lag0(ytemp,j)];
    end
    xtemp=[xtemp ones(rows(xtemp),1)];
    ytemp=ytemp(p+1:end,:);
    xtemp=xtemp(p+1:end,:);
    btemp=xtemp\ytemp;
    etemp=ytemp-xtemp*btemp;
    stemp=etemp'*etemp/(size(ytemp,1)-(p+1));
    delta(i,1) = btemp(1);
    sigma_sq(i,1)= stemp;
end
sigma = sqrt(sigma_sq);