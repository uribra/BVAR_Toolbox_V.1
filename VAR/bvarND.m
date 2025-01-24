% Written by:
% Uriel Braham
% Univerity of Bayreuth
% Universit‰tsstraﬂe 30
% 95445 Bayreuth
% uriel.braham@googlemail.com

function [resultND, VAR]    = bvarND(VARoption,BVARoption,VAR,DATA);

BVARprior.prior             = 03;
% Posterior Simulation for Coefficients with Gibbs sampling
[resultND, VAR]            = bvarNDGibbs(VARoption,BVARoption,VAR,DATA);
% Estimate the model:
[resultND, VAR]            = bvarNDpost(VARoption,BVARoption,resultND,VAR,DATA);
% Posterior Simulation for IRFs
[irf_recordND, VAR]        = bvarIRFsimul(VARoption,VAR,resultND);
% Posterior Simulation for FEVD
[fevd_recordND, VAR]       = bvarFEVDsimul(VARoption,VAR,resultND);

% Display estimation results:
bvardisplay(BVARoption,BVARprior,VARoption,VAR,resultND)
% Plot Time series
plotTS(VARoption,DATA);
% Plot  fitted values
plotFIT(VARoption,DATA,resultND);
% Plot IRF
[resultND, VAR]        = plotIRF(VARoption,VAR,resultND,irf_recordND);
% Plot Forecast Error Variance Decomoposition
[resultND, VAR] = plotFEVD(VARoption,VAR,resultND,fevd_recordND);