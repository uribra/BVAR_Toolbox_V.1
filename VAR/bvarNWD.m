% Written by:
% Uriel Braham
% Univerity of Bayreuth
% Universit‰tsstraﬂe 30
% 95445 Bayreuth
% uriel.braham@googlemail.com

function [resultNWD, VAR]   = bvarNWD(VARoption,BVARoption,VAR,DATA);

BVARprior.prior             = 05;
% Estimate the model:
[resultNWD, VAR]             = bvarNWDpost(VARoption,BVARoption,VAR,DATA);
% Posterior Simulation for Coefficients with Gibbs sampling
[resultNWD, VAR]             =  bvarNWDGibbs(VARoption,BVARoption,VAR,resultNWD);
% Posterior Simulation for IRFs
[irf_recordNWD, VAR]         = bvarIRFsimul(VARoption,VAR,resultNWD);
% Posterior Simulation for FEVD
[fevd_recordNWD, VAR]        = bvarFEVDsimul(VARoption,VAR,resultNWD);

% Display estimation results:
bvardisplay(BVARoption,BVARprior,VARoption,VAR,resultNWD)
% Plot Time series
plotTS(VARoption,DATA);
% Plot  fitted values
plotFIT(VARoption,DATA,resultNWD);
% Plot IRF
[resultNWD, VAR]             = plotIRF(VARoption,VAR,resultNWD,irf_recordNWD);
% Plot Forecast Error Variance Decomoposition
[resultNWD, VAR]             = plotFEVD(VARoption,VAR,resultNWD,fevd_recordNWD);

