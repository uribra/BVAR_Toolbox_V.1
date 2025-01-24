% Written by:
% Uriel Braham
% uriel.braham@googlemail.com

clear all; 
clear session; 
close all;  
clc; 
warning off all; 
randn('seed',2); 

%--------------------------------------------------------------------------
% SET THE DIRECTORIES
cd '/BVAR_Toolbox_V.1';
addpath(cd,'Utilities');
addpath(cd,'Subroutines');
addpath(cd,'VAR');
addpath(cd,'Data');
%--------------------------------------------------------------------------

%**************************************************************************
% OPTION SETTING FOR VAR
%**************************************************************************
%--------------------------------------------------------------------------
% EXAMPLE DATASET is from Stock and Watson (2001, JEP)
[xlsdata,xtstext]   = xlsread('data_.xls', 'Sheet1');
DATA.VARS           = xlsdata(1:end,1:end);                                     % Variables
DATA.Series         = xtstext(1,2:end);                                         % Series names
DATA.Dates          = datetime(xtstext(2:end,1) ,'InputFormat','MM\dd\yyyy');   % Note: Toolbox requires format for dates of 'MM\dd\yyyy'.
%--------------------------------------------------------------------------
% GENERAL TIME (SAMPLE) OPTIONS (structure)
VARoption.frequency             = 'q';                              % frequency of the data 'm' monthly, 'q' quarterly, 'y' annual                                         
VARoption.starty                = 1971;                             % Start year of the sample             
VARoption.startp                = 1;                                % Start period: 1=q1, 2=q2, 3=q3, 4=q4 or 1=m1,...
VARoption.endy                  = 2000;                             % End year of sample period
VARoption.endp                  = 4;                                % End period:1=q1, 2=q2, 3=q3, 4=q4 or 1=m1,...
%--------------------------------------------------------------------------
% GENERAL VAR OPTIONS (structure)
VARoption.constant              = 1;                                % 1:constant, 0=no constant
VARoption.trend                 = 0;                                % 1=linear trend, 0= linear trend
VARoption.p                     = 4;                                % number of lags
VARoption.bandwidthCoef         = 95;                               % 68=68% credible set with 16-th and 84-th percentiles, 95=95% crdible set with 
                                                                    % 2.5-th and 97.5-th percentiles
                                                                    
%--------------------------------------------------------------------------
% GIBBS SAMPLER OPTIONS (structure)
VARoption.nsave                 = 1000;
VARoption.nburn                 = 500;

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
                                                                    
%--------------------------------------------------------------------------
% FORECAST OPTION (structure)
VARoption.trainstarty           = [];                               % Start year of the training sample
VARoption.trainstartp           = [];                               % Start period of the training sample
VARoption.forhor                = [];                               % Forecast horizon in periods

%--------------------------------------------------------------------------
% OPTIONS FOR BAYESIAN PRIORS (HYPERPARAMETERS) (structure)
BVARoption.lambda1              = 0.5;                              % Overall tightness parameter
BVARoption.lambda2              = 0.5;                              % Cross-variable weighting parameter
BVARoption.lambda3              = 1;                                % Lag-decay parameter
BVARoption.lambda4              = 10000;                            % Exogenous variable tightness parameter
BVARoption.lambda5              = 10*BVARoption.lambda1;            % Sum-of-Coefficients tightness parameter
BVARoption.lambda6              = 0.2;                              % Dummy initial observation tightness parameter
BVARoption.socpri               = 0;                                % 0=no "sum-of-coefficients" prior, 1="sum-of-coefficients" prior 
BVARoption.diobs                = 0;                                % 0=no "dummy initial observation" prior, 1="dummy initial observation" prior        
BVARoption.gridl                = 0;                                % 0=no grid search for lambda1, 1=grid search for lambda1
BVARoption.marglikl             = 0;                                % 0=do not choose lambda1 from maximizing the marginal likelihood 
BVARoption.stationary           = {'nonstationary', 'nonstationary', 'nonstationary'};  % 'stationary'=stationary variable, 'nonstationary'=non-stationary variable

%--------------------------------------------------------------------------
% ENDOGENOUS VARIABLES OPTIONS (structure)
VARoption.names_endo            = {'Inflation', 'Unemployment', 'Fed Funds'}; 
VARoption.units                 = {'rate', 'rate', 'rate'};           % 'rate'= Variable is a rate, 'level'=Variable is a log-level. Note: Toolbox can only handle these two at the moment.      
%--------------------------------------------------------------------------
% EXOGENOUS VARIABLES OPTIONS (structure)
VARoption.names_exo = {};
%--------------------------------------------------------------------------
% CALL THE "OPTION' and 'GEN' FUNCTION
VARoption       = VARopt(VARoption);
[VAR, DATA]     = VARgen(VARoption,DATA);
%--------------------------------------------------------------------------
%**************************************************************************

% LAG-SELECTION:
%**************************************************************************
%pmax = 8;
%[AIChat,SIChat,HQChat]=lagorder(VARoption,DATA,pmax);
%VARoption.p = SIChat; 
%**************************************************************************


%*************************************************************************%
%                       LIST OF BAYESIAN PRIORS                           %
%*************************************************************************%
% 0.) prior = 00: (yet to be programmed, no Bayesian Approach)
%       Standard OLS VAR
% 1.) prior = 01: 
%       BVAR with diffuse ("Jeffrey's") prior on coefficients and
%       variance-covariance matrrix.
% 2.) prior = 02: (yet to be programmed)
%        Original Minnesota prior following Litterman et al. (1984)
% 3.) prior = 03; 
%       BVAR with independent Normal-Diffuse prior with Minnesota moments
% 4.) prior = 04:
%       BVAR as Normal-inverted-Wishart prior with Minnesota moments
%       following Kadiyala and Karlsson (1997)
% 5.) prior = 05:
%        BVAR as Normal-inverted-Wishart prior with Minnesota moments 
%        implemented with "dummy obervations" following Banbura et. al (2010)     
%**************************************************************************

% Choose prior:
prior = 01;

if prior == 00
    [result, VAR]   = varOLS(VARoption,DATA);
elseif prior == 01
    [result, VAR]   = bvarDiffuse(VARoption,BVARoption,VAR,DATA);
elseif prior == 03
   [result, VAR]    = bvarND(VARoption,BVARoption,VAR,DATA);
elseif prior == 04
    [result, VAR]   = bvarNW(VARoption,BVARoption,VAR,DATA);
elseif prior ==05
    [result, VAR]   = bvarNWD(VARoption,BVARoption,VAR,DATA);
end

%**************************************************************************
% END
%**************************************************************************







