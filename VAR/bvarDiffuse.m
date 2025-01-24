% Written by:
% Uriel Braham
% Univerity of Bayreuth
% Universit‰tsstraﬂe 30
% 95445 Bayreuth
% uriel.braham@googlemail.com

function [resultDiffuse, VAR]   = bvarDiffuse(VARoption,BVARoption,VAR,DATA);

BVARprior.prior                 = 01;

% Estimate the model:
[resultDiffuse, VAR]            = bvarDiffusepost(VARoption,DATA);
% Posterior Simulation for Coefficients with Gibbs sampling
[resultDiffuse, VAR]            = bvarDiffuseGibbs(VARoption,VAR,resultDiffuse);
% Posterior Simulation for IRFs
[irf_recordDiffuse, VAR]        = bvarIRFsimul(VARoption,VAR,resultDiffuse);
% Posterior Simulation for FEVD
[fevd_recordDiffuse, VAR]       = bvarFEVDsimul(VARoption,VAR,resultDiffuse);

% Display estimation results:
bvardisplay(BVARoption,BVARprior,VARoption,VAR,resultDiffuse)
% Plot Time series
plotTS(VARoption,DATA);
% Plot  fitted values
plotFIT(VARoption,DATA,resultDiffuse);
% Plot IRF
[resultDiffuse, VAR]            = plotIRF(VARoption,VAR,resultDiffuse,irf_recordDiffuse);
% Plot Forecast Error Variance Decomoposition
[resultDiffuse, VAR]            = plotFEVD(VARoption,VAR,resultDiffuse,fevd_recordDiffuse);