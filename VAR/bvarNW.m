% Written by:
% Uriel Braham
% Univerity of Bayreuth
% Universit‰tsstraﬂe 30
% 95445 Bayreuth
% uriel.braham@googlemail.com

function [resultNW, VAR]    = bvarNW(VARoption,BVARoption,VAR,DATA);

BVARprior.prior             = 04;
% Estimate the model:
[resultNW, VAR]             = bvarNWpost(VARoption,BVARoption,VAR,DATA);
% Posterior Simulation for Coefficients with Gibbs sampling
[resultNW, VAR]             = bvarNWGibbs(VARoption,BVARoption,VAR,resultNW);
% Posterior Simulation for IRFs
[irf_recordNW, VAR]         = bvarIRFsimul(VARoption,VAR,resultNW);
% Posterior Simulation for FEVD
[fevd_recordNW, VAR]        = bvarFEVDsimul(VARoption,VAR,resultNW);

% Display estimation results:
bvardisplay(BVARoption,BVARprior,VARoption,VAR,resultNW)
% Plot Time series
plotTS(VARoption,DATA);
% Plot  fitted values
plotFIT(VARoption,DATA,resultNW);
% Plot IRF
[resultNW, VAR]             = plotIRF(VARoption,VAR,resultNW,irf_recordNW);
% Plot Forecast Error Variance Decomoposition
[resultNW, VAR]             = plotFEVD(VARoption,VAR,resultNW,fevd_recordNW);