function VARoption = VARopt(VARoption)

%**************************************************************************
% OPTION SETTING FOR VAR
%**************************************************************************

%--------------------------------------------------------------------------
% GENERAL VAR OPTIONS (structure)
VARoption.bandwidthCoef         = 95;                               % 68=68% credible set with 16-th and 84-th percentiles, 95=95% crdible set with 
                                                                    % 2.5-th and 97.5-th percentiles
%--------------------------------------------------------------------------
% GIBBS SAMPLER OPTIONS (structure)
VARoption.nsave                 = 20000;
VARoption.nburn                 = 5000;
%--------------------------------------------------------------------------
% IMPLUSE RESPONSE OPTIONS (structure)
VARoption.ihor                  = 24;                               % Horizon for Impluse Responses
VARoption.shocksize             = 1;                                % 0=standard deviation shock, 1=shock normalized to unity
VARoption.bandwidthIRF          = 68;                               % 68=68% credible set with 16-th and 84-th percentiles, 95=95% crdible set with 
                                                                    % 2.5-th and 97.5-th percentiles                                                                 
%--------------------------------------------------------------------------
% FORECAST ERROR VARIANCE DECOMPOISITION OPTIONS (structure)  
VARoption.ihorfevd              = 24;                               % Horizon for Forecast Error Variance Decomposition
VARoption.bandwidthFEVD         = 68;                               % 68=68% credible set with 16-th and 84-th percentiles, 95=95% crdible set with 
                                                                    % 2.5-th and 97.5-th percentiles