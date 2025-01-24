% Written by:
% Uriel Braham
% uriel.braham@googlemail.com

function [lambdaopt] = bvarNWDgridsearch(VARoption,VAR,DATA)

GRID = [0:.025:5 50];

GRIDsearch = GRID*sqrt(VAR.M*VARoption.p);
logmarglik = zeros(1,length(GRIDsearch));
for jpi = 1:1:length(GRIDsearch)
    
    VARoption.lambda1 = 1/GRIDsearch(jpi);   % The overall tightness prior
    VARoption.lambda5 = 10*VARoption.lambda1;             % Set the prior on the sum of the coefficients
    % Estimate the Bayesian VAR to get the fit corresponding to the prior pi
    [result, ~] = bvarNWD(VARoption,DATA);
    logmarglik(1,jpi) = result.lnmarglik;
end

% Find lambda maximizing log marginal likelihood
[temp,Jstar]    = max(-(logmarglik));
logmarglikstar  = logmarglik(Jstar);
il1             = GRIDsearch(Jstar);
lambdaopt       = 1/il1;