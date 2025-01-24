% Written by:
% Uriel Braham
% uriel.braham@googlemail.com

function bvardisplay(BVARoption,BVARprior,VARoption,VAR,result)

%**************************************************************************
%   FUNCTION TO DISPLAY THE RESULTS OF AN ESTIMATED BAYESIAN VAR MODEL
%--------------------------------------------------------------------------
%   INPUT: 
%   - Structure 'BVARoption'
%   - Structure 'BVARprior'
%   - Structure 'VARoption'
%   - Structure 'VARoption'
%   - Structure 'VAR'
%   - Structure 'result'
%--------------------------------------------------------------------------
%   OUTPUT: 
%   - Display of estimation information in the console and .txt file
%   - Display of estimation results in the console
%--------------------------------------------------------------------------
%
%**************************************************************************

%----------------------------CHECK INPUTS----------------------------------
if ~exist('BVARoption')
    error('You need to provide BVAR options in "BVARoption".');
end
if ~exist('VARoption')
    error('You need to provide VAR options in "VARoption".');
end
if ~exist('BVARprior')
    error('You need to provide the prior for the BVAR in "BVARprior".');
end
if ~exist('VAR')
    error('You need to provide estmation information in "VAR".');
end
if ~exist('result')
    error('You need to provide estimation results in "results".');
end

%--------------------------------------------------------------------------
% Define the ".txt" file to save the estimation results
filelocation=[cd '/Results/estimationresults.txt'];
savefile=fopen(filelocation,'w+');
%--------------------------------------------------------------------------

%**************************************************************************
% UNPACK structures 'BVARoption','BVARprior','VARoption' and 'VAR'
%**************************************************************************

lambda1         = BVARoption.lambda1;
lambda2         = BVARoption.lambda2;
lambda3         = BVARoption.lambda3;
lambda4         = BVARoption.lambda4;
lambda5         = BVARoption.lambda5;
lambda6         = BVARoption.lambda6;
socpri          = BVARoption.socpri; 
diobs           = BVARoption.diobs;
stationary      = BVARoption.stationary;
prior           = BVARprior.prior;
%
frequency       = VARoption.frequency;                                                                               
starty          = VARoption.starty;                                  
startp          = VARoption.startp;                                        
endy            = VARoption.endy;                                           
endp            = VARoption.endp;
datesarray      = DatesOption(starty,endy,frequency,startp,endp);
stringdates     = datesarray.dates;
names_endo      = VARoption.names_endo;
names_exo       = VARoption.names_exo;
constant        = VARoption.constant;
trend           = VARoption.trend;
p               = VARoption.p;
bandwidthCoef   = VARoption.bandwidthCoef; 
bandwidthIRF    = VARoption.bandwidthIRF;

nsave           = VARoption.nsave;   
nburn           = VARoption.nburn;         
ihor            = VARoption.ihor;  
shocksize       = VARoption.shocksize; 

K               = VAR.K;
M               = VAR.M;
k               = VAR.k;
N               = VAR.N;
q               = VAR.q;
T               = VAR.T;

%**************************************************************************
% UNPACK structure 'result'
%**************************************************************************

A_post          = result.A_post;
a_post          = result.a_post;
%a_SD           = result.a_SD;
A_draws         = result.A_draws;
a_draws         = result.a_draws;
SIGMA_post      = result.SIGMA_post;
SIGMA_draws     = result.SIGMA_draws;
SSE             = result.SSE;
RSME            = result.RMSE;
R2              = result.R2;
S               = result.S;
maxeig          = result.maxeig;
%FIT             = result.FIT;

%*************************************************************************%
%                          ESTIMATION RESULTS                             %
%*************************************************************************%
clc;
fprintf('%s\n', 'Written by: '); 
fprintf(savefile,'%s\n', 'Written by: ');
fprintf('%s\n', 'Uriel Braham'); 
fprintf(savefile,'%s\n', 'Uriel Braham');
fprintf('%s\n', 'uriel.braham@googlemail.com');
fprintf(savefile,'%s\n','uriel.braham@googlemail.com');
fprintf('%s\n', '');
fprintf(savefile,'%s\n','');

fprintf('%s\n', 'Disclaimer');
fprintf(savefile,'%s\n','Disclaimer');
fprintf('%s\n', 'The code of this toolbox can be run and rewritten without my permission.');
fprintf(savefile,'%s\n','The code of this toolbox can be run and rewritten without my permission.');
fprintf('%s\n', 'I will not be responsible for any damage/liability.');
fprintf(savefile,'%s\n','I will not be responsible for any damage/liability.');
fprintf('%s\n', '');
fprintf(savefile,'%s\n','');

time=clock;
datestring=datestr(time);
dateinfo=['Date: ' datestring(1,1:11) '   Time: ' datestring(1,13:17)];
fprintf('%s\n',dateinfo);
fprintf(savefile,'%s\n',dateinfo);
fprintf(savefile,'%s\n','');

fprintf('%s\n','%---------------------------------------------------------------------------------------%');
fprintf(savefile,'%s\n','%---------------------------------------------------------------------------------------%');
fprintf('%s\n','%                                ESTIMATION RESULTS                                     %');
fprintf(savefile,'%s\n','%                                ESTIMATION RESULTS                                     %');
fprintf('%s\n','%---------------------------------------------------------------------------------------%');
fprintf(savefile,'%s\n','%---------------------------------------------------------------------------------------%');

%*************************************************************************%
%                          BVAR INFORMATION                               %
%*************************************************************************%
fprintf('%s\n','');
fprintf(savefile,'%s\n','');
fprintf('%s\n','Bayesian VAR information ');
fprintf(savefile,'%s\n','Bayesian VAR information');
fprintf('%s\n','-----------------------------------------------------------------------------------------');
fprintf(savefile,'%s\n','-----------------------------------------------------------------------------------------');
% Endogenous variables
temp='Endogenous variables: ';
for ii=1:M
temp=[temp ' ' names_endo{1,ii} ' '];
end
endoinfo=temp;
fprintf('%s\n',['Number of endogenous variables: ' num2str(M)]);
fprintf('%s\n',endoinfo);
fprintf(savefile,'%s\n',endoinfo);
% Exogenous variables
temp='Exogenous variables: ';
if constant==0 && trend== 0 && q==0
   temp=[temp ' None'];
elseif constant==1 && trend == 1 && q==2
   temp=[temp ' Constant ' ' Linear trend '];
elseif constant==1 && trend == 0 && q==1
    temp = [temp ' Constant'];
elseif constant==0 && trend ==1 && q==1
    temp = [temp ' Linear trend'];
elseif constant==0 && trend == 0 &&  q>0
   temp=[];
   for ii=1:N
        temp=[temp ' ' names_exo{1,ii} ' '];
   end
elseif constant==0 && trend == 1 &&  q>1
   temp=[temp ' Linear trend'];
   for ii=1:N
        temp=[temp ' ' names_exo{1,ii} ' '];
   end  
elseif constant==1 && trend == 1 &&  q>2
   temp=[temp ' Constant' ' Linear trend'];
   for ii=1:N
        temp=[temp ' ' names_exo{1,ii} ' '];
   end     
elseif constant==1 && trend == 0 &&  q>1  
   temp=[temp ' Constant'];
   for ii=1:N
        temp=[temp ' ' names_exo{1,ii} ' '];
   end     
end
exoinfo=temp;
fprintf('%s\n',['Number of exogenous variables: ' num2str(q)]);
fprintf('%s\n',exoinfo);
fprintf(savefile,'%s\n',exoinfo);

if frequency == 'q'
    freq_label = 'Quarters';
elseif frequency == 'm'
    freq_label = 'Months';
elseif VARoption.frequency == 'y'
    freq_label = 'Years';
end

%Frequency
freq = ['Data frequency: ' freq_label];
fprintf('%s\n',freq);
fprintf(savefile,'%s\n',freq);

% SAMPLE PERIOD 
startdate = stringdates{1};
enddate   = stringdates{end};
sampledateinfo=['Estimation sample period: ' startdate '-' enddate];
fprintf('%s\n',sampledateinfo);
fprintf(savefile,'%s\n',sampledateinfo);
% Sample size 
samplelengthinfo=['Sample size (omitting initial conditions): ' num2str(T)];
fprintf('%s\n',samplelengthinfo);
fprintf(savefile,'%s\n',samplelengthinfo);
% Number of lag order on 
laginfo=['Number of autoregressive lags: ' num2str(p)];
fprintf('%s\n',laginfo);
fprintf(savefile,'%s\n',laginfo);

fprintf('%s\n','');


% PRIOR CHOSEN
if prior==01
priorused='Prior: Jeffreys Prior';
elseif prior == 02
priorused='Prior: Minnesota (SIGMA as full VAR estimates)';
elseif prior==03
priorused='Prior: Normal-Diffuse';
elseif prior==04
priorused='Prior: Normal-Wishart';
elseif prior==05
priorused='Prior: Dummy Observations';
end
fprintf('%s\n',priorused);
fprintf(savefile,'%s\n',priorused);

% HYPERPARAMETERS
hyperparam = 'Hyperparameters:';
fprintf('%s\n',hyperparam);
fprintf(savefile,'%s\n',hyperparam);

if prior==01
    fprintf('%s\n','No hyperparameter to select.');
    fprintf(savefile,'%s\n','No hyperparameter to select.');
end

delta    = ones(length(stationary),1); 
delta(ismember(stationary, 'stationary')) = 0;  
deltai              = delta'; 
arc = cell(1,M);
for iii=1:M
    if deltai(1,iii) == 1
       arc{1,iii} = 'Random Walk';
    elseif deltai(1,iii) == 0
        arc{1,iii}  = 'White Noise';
    end
end 
temp='Autoregressive Coefficient: ';
for ii=1:3
temp=[temp ' ' arc{1,ii} ' '];
end
if prior==02 || prior ==03 || prior ==04 || prior ==05
hyperparameter2=temp;
fprintf('%s\n',hyperparameter2);
fprintf(savefile,'%s\n',hyperparameter2);
end

if prior==02 || prior ==03 || prior ==04 || prior ==05
hyperparameter3=['Overall tightness (lambda1):                    ' num2str(lambda1)];
fprintf('%s\n',hyperparameter3);
fprintf(savefile,'%s\n',hyperparameter3);
end

if prior==02
hyperparameter4=['Cross-Variable weighting (lambda2):             ' num2str(lambda2)];
fprintf('%s\n',hyperparameter4);
fprintf(savefile,'%s\n',hyperparameter4);
end

if prior==02 || prior ==03 || prior ==04 || prior ==05
hyperparameter5=['Lag decay (lambda3):                            ' num2str(lambda3)];
fprintf('%s\n',hyperparameter5);
fprintf(savefile,'%s\n',hyperparameter5);
end

if prior==02 || prior ==03 || prior ==04 || prior ==05
hyperparameter6=['Exogenous variable tightness (lambda4):         ' num2str(lambda4)];
fprintf('%s\n',hyperparameter6);
fprintf(savefile,'%s\n',hyperparameter6);
end

if  (prior==02 || prior ==03 || prior ==04 || prior ==05) && socpri==1
hyperparameter7=['Sum-of-Coefficients tightness (lambda5):        ' num2str(lambda5)];
fprintf('%s\n',hyperparameter7);
fprintf(savefile,'%s\n',hyperparameter7);
end

if (prior==02 || prior ==03 || prior ==04 || prior ==05) && diobs==1
hyperparameter8=['Dummy-initial-observation tightness (lambda6):  ' num2str(lambda6)];
fprintf('%s\n',hyperparameter8);
fprintf(savefile,'%s\n',hyperparameter8);
end

% GIBBS SAMPLER INFORMATION
gibbsinfo = 'Gibbs sampler information: ';
fprintf('%s\n',gibbsinfo);
fprintf(savefile,'%s\n',gibbsinfo);
gibbssave=['Number of saved Gibbs draws: ' num2str(nsave)];
fprintf('%s\n',gibbssave);
fprintf(savefile,'%s\n',gibbssave);
gibbsburn=['Number of burn-in Gibbs draws: ' num2str(nburn)];
fprintf('%s\n',gibbsburn);
fprintf(savefile,'%s\n',gibbsburn);

% IMPLUSE RESPONSE INFORMATION
IRFinfo = 'Impulse Response Function (IRF) information: ';
fprintf('%s\n',IRFinfo);
fprintf(savefile,'%s\n',IRFinfo);
IRFihor=['IRF periods: ' num2str(ihor)];
fprintf('%s\n',IRFihor);
fprintf(savefile,'%s\n',IRFihor);

if shocksize             == 1
        IRFshock = 'Triangular factorization (Cholesky decomposition)';
elseif shocksize         == 0
        IRFshock = 'One standard deviation (Cholesky decomposition)';
end
fprintf('%s\n',['Structural shock size: ' IRFshock]);
fprintf(savefile,'%s\n',['Structural shock size: ' IRFshock]);
IRFband=['IRF credible set (HPDI): ' num2str(bandwidthIRF) ' Percent'];
fprintf('%s\n',IRFband);
fprintf(savefile,'%s\n',IRFband);
fprintf('%s\n','-----------------------------------------------------------------------------------------');
fprintf(savefile,'%s\n','-----------------------------------------------------------------------------------------');
fprintf('%s\n',' ');
fprintf(savefile,'%s\n','');
fprintf('%s\n',' ');
fprintf(savefile,'%s\n','');

% STABILITY OF THE VAR
fprintf('%s\n','Stability ');
fprintf(savefile,'%s\n','Stability ');
fprintf('%s\n','-----------------------------------------------------------------------------------------');
fprintf(savefile,'%s\n','-----------------------------------------------------------------------------------------');
stabilityinfo1=['Maximum Eigenvalue of the companion matrix: ' num2str(maxeig)];
fprintf('%s\n',[stabilityinfo1]);
fprintf(savefile,'%s\n',stabilityinfo1);
if S==1;
stabilityinfo2=['All Eigenvalues lie within the unit circle.'];
stabilityinfo21=['The estimated VAR model satisfies the stability condition.'];
fprintf('%s\n',stabilityinfo2);
fprintf(savefile,'%s\n',stabilityinfo2);
fprintf('%s\n',stabilityinfo21);
fprintf(savefile,'%s\n',stabilityinfo21);
else
stabilityinfo2=['Warning: At least one Eigenvalue lies outside the unit circle.'];
stabilityinfo21=['The estimated VAR model is not be stable.'];
fprintf('%s\n',stabilityinfo2);
fprintf(savefile,'%s\n',stabilityinfo2);
fprintf('%s\n',stabilityinfo21);
fprintf(savefile,'%s\n',stabilityinfo21);
end
fprintf('%s\n','-----------------------------------------------------------------------------------------');
fprintf(savefile,'%s\n','-----------------------------------------------------------------------------------------');


%*************************************************************************%
%                           BVAR ESTIMATES                                %
%*************************************************************************%
fprintf('%s\n','');
fprintf(savefile,'%s\n','');
fprintf('%s\n','');
fprintf(savefile,'%s\n','');

lower_q = ((100-bandwidthCoef)/2)/100; 
upper_q = (100 - (100-bandwidthCoef)/2)/100;

lower = reshape(quantile(a_draws,lower_q)',K,M);
upper = reshape(quantile(a_draws,upper_q)',K,M);
A_std = reshape(std(a_draws)',K,M);

coeffinfo=['Reduced-form BVAR coefficients: Posterior estimates'];
fprintf('%s\n',coeffinfo);
fprintf(savefile,'%s\n',coeffinfo);
fprintf('%s\n','-----------------------------------------------------------------------------------------');
fprintf(savefile,'%s\n','-----------------------------------------------------------------------------------------');
bandinfo = ['Coefficient credible set (HPDI): ' num2str(bandwidthCoef) ' Percent'];
fprintf('%s\n',bandinfo);
fprintf(savefile,'%s\n',bandinfo);


for ii=1:M % loop over endogenous variables 

    fprintf('%s\n','');
    fprintf(savefile,'%s\n','');
    if ii~=1
        fprintf('%s\n','');
    fprintf(savefile,'%s\n','');
    end

    endoinfo=['Endogenous variable: ' names_endo{1,ii}];
    fprintf('%s\n',endoinfo);
    fprintf(savefile,'%s\n',endoinfo);
    fprintf('%s\n','-----------------------------------------------------------------------------------------');
    fprintf(savefile,'%s\n','-----------------------------------------------------------------------------------------');
    header=fprintf('%25s %15s %15s %15s %15s\n','','Median','St.Dev','Lower bound','Upper bound');
    header=fprintf(savefile,'%25s %15s %15s %15s %15s\n','','Median','St.Dev','Lower bound','Upper bound');
    fprintf('%s\n','-----------------------------------------------------------------------------------------');
    fprintf(savefile,'%s\n','-----------------------------------------------------------------------------------------');

   ccc = 0; 
   for jj=1:p % loop over lags 
       for kk=1:M % loop over endogenous variables
           ccc= ccc+1;
           valuesendo = [A_post(ccc,ii) A_std(ccc,ii) lower(ccc,ii) upper(ccc,ii)];
           fprintf('%25s %15.3f %15.3f %15.3f %15.3f\n',strcat(names_endo{1,kk},'(t-',int2str(jj),')'),valuesendo);
           fprintf(savefile,'%25s %15.3f %15.3f %15.3f %15.3f\n',strcat(names_endo{1,kk},'(t-',int2str(jj),')'),valuesendo);
       end 
   end 
   
   if constant == 1 && trend ==1 
       valuesexo1 = [A_post(M*p+1:M*p+2,ii) A_std(M*p+1:M*p+2,ii) lower(M*p+1:M*p+2,ii) upper(M*p+1:M*p+2,ii)];
       fprintf('%25s %15.3f %15.3f %15.3f %15.3f\n',strcat(['Constant'; 'Linear Trend']),valuesexo1);
       fprintf(savefile,'%25s %15.3f %15.3f %15.3f %15.3f\n',strcat(['Constant'; 'Linear Trend']),valuesexo1);
       if N > 0
           for ll=1:N %loop over exogenous variables
                valuesexo2 = [A_post(M*p + constant + trend+1:end ,ii) A_std(M*p + constant + trend:end ,ii) lower(M*p + constant + trend:end ,ii) upper(M*p + constant + trend:end ,ii)];
                fprintf('%25s %15.3f %15.3f %15.3f %15.3f\n',strcat(names_exo{1,ll}),valuesexo2);
                fprintf(savefile,'%25s %15.3f %15.3f %15.3f %15.3f\n',strcat(names_exo{1,ll}),valuesexo2);
           end  % end loop over exogenous variables
       end 
       
   elseif constant == 1 && trend ==0
       valuesexo1 = [A_post(M*p+1,ii) A_std(M*p+1,ii) lower(M*p+1,ii) upper(M*p+1,ii)];
       fprintf('%25s %15.3f %15.3f %15.3f %15.3f\n',strcat(['Constant']),valuesexo1);
       fprintf(savefile,'%25s %15.3f %15.3f %15.3f %15.3f\n',strcat(['Constant']),valuesexo1);
       if N > 0
           for ll=1:N %loop over exogenous variables
                valuesexo2 = [A_post(M*p + constant + trend+1:end ,ii) A_std(M*p + constant + trend:end ,ii) lower(M*p + constant + trend:end ,ii) upper(M*p + constant + trend:end ,ii)];
                fprintf('%25s %15.3f %15.3f %15.3f %15.3f\n',strcat(names_exo{1,ll}),valuesexo2);
                fprintf(savefile,'%25s %15.3f %15.3f %15.3f %15.3f\n',strcat(names_exo{1,ll}),valuesexo2);
           end  % end loop over exogenous variables
       end 
       
   elseif constant == 0 && trend ==1
       valuesexo1 = [A_post(M*p+1,ii) A_std(M*p+1,ii) lower(M*p+1,ii) upper(M*p+1,ii)];
       fprintf('%25s %15.3f %15.3f %15.3f %15.3f\n',strcat(['Linear Trend']),valuesexo1);
       fprintf(savefile,'%25s %15.3f %15.3f %15.3f %15.3f\n',strcat(['Linear Trend']),valuesexo1);
       if N > 0
           for ll=1:N %loop over exogenous variables
                valuesexo2 = [A_post(M*p + constant + trend+1:end ,ii) A_std(M*p + constant + trend:end ,ii) lower(M*p + constant + trend:end ,ii) upper(M*p + constant + trend:end ,ii)];
                fprintf('%25s %15.3f %15.3f %15.3f %15.3f\n',strcat(names_exo{1,ll}),valuesexo2);
                fprintf(savefile,'%25s %15.3f %15.3f %15.3f %15.3f\n',strcat(names_exo{1,ll}),valuesexo2);
           end  % end loop over exogenous variables
       end 
        
   elseif constant == 0 && trend ==0
       
       if N > 0
           for ll=1:N %loop over exogenous variables
                valuesexo2 = [A_post(M*p + constant + trend+1:end ,ii) A_std(M*p + constant + trend:end,ii) lower(M*p + constant + trend:end ,ii) upper(M*p + constant + trend:end ,ii)];
                fprintf('%25s %15.3f %15.3f %15.3f %15.3f\n',strcat(names_exo{1,ll}),valuesexo2);
                fprintf(savefile,'%25s %15.3f %15.3f %15.3f %15.3f\n',strcat(names_exo{1,ll}),valuesexo2);
           end   % end loop over exogenous variables
       end   
   end
    
   fprintf('%s\n','');
   fprintf(savefile,'%s\n','');

%*************************************************************************%
%                       FURTHER STATISTICS                                %
%*************************************************************************%

   sse = ['Sum of squared residuals: ' num2str(SSE(ii,ii),'%.2f')];
   fprintf('%s\n',sse);
   fprintf(savefile,'%s\n',sse);

   rsquared =['R-squared: ' num2str(R2(ii,1),'%.3f')];
   fprintf('%s\n',rsquared);
   fprintf(savefile,'%s\n',rsquared);
    
   rmse = ['Root Mean Squared Error (RMSE): ' num2str(RSME(ii,1),'%.3f')];
   fprintf('%s\n',rmse);
   fprintf(savefile,'%s\n',rmse);
   
fprintf('%s\n','-----------------------------------------------------------------------------------------');
fprintf(savefile,'%s\n','-----------------------------------------------------------------------------------------');
end

fprintf('%s\n','');
fprintf(savefile,'%s\n','');

SIGMAinfo=['SIGMA (residual Variance-Covariance matrix): Posterior estimates'];
fprintf('%s\n',SIGMAinfo);
fprintf(savefile,'%s\n',SIGMAinfo);
width=length(sprintf('%d',floor(max(abs(vec(SIGMA_post))))));
width=width+5;
for ii=1:M
tempSIG=[];
   for jj=1:M
   % convert matrix entry into string
   number=num2str(SIGMA_post(ii,jj),'% .3f');
      % pad potential missing blanks
      while numel(number)<width
      number=[' ' number];
      end
   number=[number '  '];
   tempSIG=[tempSIG number];
   end
fprintf('%s\n',tempSIG);
fprintf(savefile,'%s\n',tempSIG);
end
fprintf('%s\n','-----------------------------------------------------------------------------------------');
fprintf(savefile,'%s\n','-----------------------------------------------------------------------------------------');
fprintf('%s\n','');
fprintf(savefile,'%s\n','');


%*************************************************************************%
%              POSTERIOR DISTRIBUTIONS OF COEFFICIENTS                    %
%*************************************************************************%
% Figure properties
bins = 40;
titleFontSize = 11;

for ii=1:M % loop over endogenous variables

    figname=['Posterior distribution: ' ' ' names_endo{1,ii} ' '];
    set(figure,'name', figname);
    set(gca,'Color',[0.9 0.9 0.9]);
    counterplot = 0;
    
    % Plot of posterior distribution of autoregressive coefficients
    for jj=1:p % loop over lags of endogenous variables 
        for kk=1:M % loop over endogenous variables
            counterplot = counterplot + 1; 
            subplot(p+ceil(q/M),M,counterplot) 
            valuesendo = [A_draws(:,counterplot,ii)];
            histogram(valuesendo, bins, 'Normalization','probability')
            properties=findobj(gca,'Type','patch');
            set(properties,'FaceColor',[0.7 0.78 1],'EdgeColor',[0 0 0],'LineWidth',0.1);
            title([names_endo{1,kk} '(t-' int2str(jj) ')'], 'FontSize', titleFontSize, 'FontWeight','normal');
        end   
    end  
    
    
     if constant == 1 && trend ==1
            counterplot = counterplot+1;
            valuesexo1 = [A_draws(:,M*p+1,ii)];
            subplot(p+ceil(q/M),M,counterplot) 
            histogram(valuesexo1, bins, 'Normalization','probability')
            properties=findobj(gca,'Type','patch');
            set(properties,'FaceColor',[0.7 0.78 1],'EdgeColor',[0 0 0],'LineWidth',0.1);
            title(['Constant'], 'FontSize', titleFontSize, 'FontWeight','normal');
            
            counterplot = counterplot+1;
            valuesexo1 = [A_draws(:,M*p+2,ii)];
            subplot(p+ceil(q/M),M,counterplot) 
            histogram(valuesexo1, bins, 'Normalization','probability')
            properties=findobj(gca,'Type','patch');
            set(properties,'FaceColor',[0.7 0.78 1],'EdgeColor',[0 0 0],'LineWidth',0.1);
            title(['Linear Trend'], 'FontSize', titleFontSize, 'FontWeight','normal');
            
        if N > 0
            for ll=1:N %loop over exogenous variables
                counterplot = counterplot + 1; 
                subplot(p+ceil(q/M),M,counterplot) 
                valuesexo2 = [A_post(:,M*p + constant + trend+1:end ,ii)];
                histogram(valuesexo2, bins, 'Normalization','probability')
                properties=findobj(gca,'Type','patch');
                set(properties,'FaceColor',[0.7 0.78 1],'EdgeColor',[0 0 0],'LineWidth',0.1);
                title([names_exo{1,ll}], 'FontSize', titleFontSize, 'FontWeight','normal');   
           end  % end loop over exogenous variables
        end 
     elseif constant == 1 && trend ==0 
                counterplot = counterplot+1;
                valuesexo1 = [A_draws(:,M*p+1,ii)];
                subplot(p+ceil(q/M),M,counterplot) 
                histogram(valuesexo1, bins, 'Normalization','probability')
                properties=findobj(gca,'Type','patch');
                set(properties,'FaceColor',[0.7 0.78 1],'EdgeColor',[0 0 0],'LineWidth',0.1);
                title(['Constant'], 'FontSize', titleFontSize, 'FontWeight','normal');  
            
            if N > 0
                for ll=1:N %loop over exogenous variables
                    counterplot = counterplot + 1; 
                    subplot(p+ceil(q/M),M,counterplot) 
                    valuesexo2 = [A_post(:,M*p + constant + trend+1:end ,ii)];
                    histogram(valuesexo2, bins, 'Normalization','probability')
                    properties=findobj(gca,'Type','patch');
                    set(properties,'FaceColor',[0.7 0.78 1],'EdgeColor',[0 0 0],'LineWidth',0.1);
                    title([names_exo{1,ll}], 'FontSize', titleFontSize, 'FontWeight','normal');   
                end  % end loop over exogenous variables
            end 

      elseif constant == 0 && trend ==1
            counterplot = counterplot+1;
            valuesexo1 = [A_draws(:,M*p+1,ii)];
            subplot(p+ceil(q/M),M,counterplot) 
            histogram(valuesexo1, bins, 'Normalization','probability')
            properties=findobj(gca,'Type','patch');
            set(properties,'FaceColor',[0.7 0.78 1],'EdgeColor',[0 0 0],'LineWidth',0.1);
            title(['Linear Trend'], 'FontSize', titleFontSize, 'FontWeight','normal');   
                if N > 0
                    for ll=1:N %loop over exogenous variables
                        counterplot = counterplot + 1; 
                        subplot(p+ceil(q/M),M,counterplot) 
                        valuesexo2 = [A_post(:,M*p + constant + trend+1:end ,ii)];
                        histogram(valuesexo2, bins, 'Normalization','probability')
                        properties=findobj(gca,'Type','patch');
                        set(properties,'FaceColor',[0.7 0.78 1],'EdgeColor',[0 0 0],'LineWidth',0.1);
                        title([names_exo{1,ll}], 'FontSize', titleFontSize, 'FontWeight','normal');   
                    end  % end loop over exogenous variables 
                end 
         
         
       elseif constant == 0 && trend ==0 
            
                if N > 0
                      for ll=1:N %loop over exogenous variables
                         counterplot = counterplot + 1; 
                            subplot(p+ceil(q/M),M,counterplot) 
                            valuesexo2 = [A_post(:,M*p + constant + trend+1:end ,ii)];
                            histogram(valuesexo2, bins, 'Normalization','probability')
                            properties=findobj(gca,'Type','patch');
                            set(properties,'FaceColor',[0.7 0.78 1],'EdgeColor',[0 0 0],'LineWidth',0.1);
                            title([names_exo{1,ll}], 'FontSize', titleFontSize, 'FontWeight','normal');   
                     end  % end loop over exogenous variables
                end 
         
     end
     
    ax=axes('Units','Normal','Position',[.11 .075 .85 .88],'Visible','off');
    set(get(ax,'Title'),'Visible','on')
    title(['Endogenous variable: ' names_endo{1,ii} ],'FontSize',12,...%'FontName','Times New Roman',
    'FontWeight','normal');


    %filename = ['PosteriorDist_'];
    %FigName = [filename num2str(ii)];
	%set(gcf, 'Color', 'w');
    %export_fig(FigName,'-pdf','-png','-painters')
             
    filename = ['PosteriorDist_'];
    FigName = [filename num2str(ii)];
    path = [cd '/Results'];
    saveas(gcf, fullfile(path, FigName), 'pdf')
    
end % end loo over endogenous variables



