% Written by:
% Uriel Braham
% uriel.braham@googlemail.comm

function [NDprior] = bvarNDprior(VARoption,BVARoption,DATA)

%**************************************************************************
%   FUNCTION TO DRAW FROM POSTERIOR DISTIBUTIONS WITH GIBBS SAMPLER
%--------------------------------------------------------------------------
%   INPUT: Structure 'DATA' with
%   -
%--------------------------------------------------------------------------
%   OUTPUT: STRUCTURE OF NWM RESULTS
%   -   
%--------------------------------------------------------------------------
%
%**************************************************************************

%==========================================================================
% VAR OPTIONS
%==========================================================================
constant    = VARoption.constant;
p           = VARoption.p;
trend       = VARoption.trend;

%==========================================================================
% DATA 
%==========================================================================
Yendo       = DATA.Yendo;
Yexo        = DATA.Yexo;
M           = cols(Yendo);
N           = cols(Yexo);
q           = constant + trend + N;
K           = M*p + q;
k           = M*K;

%==========================================================================
% PRIOR OPTIONS
%==========================================================================
lambda1             = BVARoption.lambda1;
lambda2             = BVARoption.lambda2;
lambda3             = BVARoption.lambda3;
lambda4             = BVARoption.lambda4;
lambda5             = BVARoption.lambda5; 
lambda6             = BVARoption.lambda6;
BVARoption.delta    = ones(length(BVARoption.stationary),1); 
BVARoption.delta(ismember(BVARoption.stationary, 'stationary')) = 0;  
deltai              = BVARoption.delta'; 

socpri              = BVARoption.socpri; 
diobs               = BVARoption.diobs;


il1                 = 1/lambda1;
il2                 = 1/lambda2;
il3                 = 1/lambda3;
il4                 = 1/lambda4;
il5                 = 1/lambda5;
il6                 = 1/lambda6;

%==========================================================================
% DATA PREPS
%==========================================================================
% Endogenous and exogenous variables
Yendo       = DATA.Yendo;
Yexo        = DATA.Yexo;
% Full sample mean:
muYendo     = mean(Yendo);
% Pre sample mean
muYendo0    = mean(Yendo(1:p,:),1);
% Estimated RMSE for each variable of AR(p)
sigma_sq    = zeros(M,1);               
pm = p;
for i = 1:M
    ytemp=Yendo(:,i);
    xtemp=[];
    for j=1:pm
    xtemp=[xtemp lag0(ytemp,j)];
    end
    xtemp=[xtemp ones(rows(xtemp),1)];
    ytemp=ytemp(pm+1:end,:);
    xtemp=xtemp(pm+1:end,:);
    btemp=xtemp\ytemp;
    etemp=ytemp-xtemp*btemp;
    stemp=etemp'*etemp/(size(ytemp,1)-(pm+1));
    sigma_sq(i,1)= stemp;
end
sigma = sqrt(sigma_sq)';


%==========================================================================
% PRIOR
%==========================================================================
A_prior         = zeros(K,M);
A_prior(1:M,:)  = eye(M)*diag(deltai);
a_prior         = vec(A_prior);

% next compute phi0, the variance-covariance matrix of beta, defined in (1.4.7)

v_prior = 0;

% set first phi0 as a k*k matrix of zeros
PHI_prior=zeros(k,k);

% set variance for coefficients on own lags, using (1.3.5)
for ii=1:M
   for jj=1:p
   PHI_prior((ii-1)*K+(jj-1)*M+ii,(ii-1)*K+(jj-1)*M+ii)=(lambda1/jj^lambda3)^2;
   end
end


%  set variance for coefficients on cross lags, using (1.3.6)
for ii=1:M
   for jj=1:p
      for kk=1:M
      if kk==ii
      else
      PHI_prior((ii-1)*K+(jj-1)*M+kk,(ii-1)*K+(jj-1)*M+kk)=(sigma_sq(ii,1)/sigma_sq(kk,1))*(((lambda1*lambda2)/(jj^lambda3))^2);
      end
      end
   end
end

% finally set the variance for exogenous variables, using (1.3.7)
for ii=1:M 
   for jj=1:q
   PHI_prior(ii*K-q+jj,ii*K-q+jj)=sigma_sq(ii,1)*((lambda1*lambda4)^2);
   end
end

%==========================================================================
% DUMMIES
%==========================================================================
% Construct dummy for the sum of the coefficients
Yd1 = il5*diag(muYendo0.*deltai);
Xd1 = [kron(ones(1,p),Yd1) zeros(M,q)];

% Construct dummy for the dummy initial observation   
% Generate the regressor matrix
if constant ==1
    if isempty(Yexo)==1
        if trend==1
                X1 = [ones(size(Yendo,1),1) (1:size(Yendo,1))'];
        elseif trend==0 
                X1 = [ones(size(Yendo,1),1)];  
         end 
    elseif isempty(Yexo)==0
         if trend ==1
                X1 = [ones(size(Yendo,1),1) (1:size(Yendo,1))' Yexo];
         elseif trend==0 
                X1 = [ones(size(Yendo,1),1) Yexo]; 
         end
    end 
elseif constant==0    
    if isempty(Yexo)==1
        if trend==1
                X1 = [(1:size(Yendo,1))'];
        elseif trend==0 
                X1 = [];  
        end 
    elseif isempty(Yexo)==0
        if trend ==1
                X1 = [(1:size(Yendo,1))' Yexo];
        elseif trend==0 
                X1 = [Yexo]; 
        end
    end  
end

if isempty(X1)
   Xbar=[];
else
   Xbar=mean(X1(1:p,:));
end
% Construct for initial observation dummy
Yd2= il6*muYendo0;  
Xd2=[kron(ones(1,p),Yd2) il6*Xbar];

% Generate dummy obervations
if BVARoption.socpri == 0 && BVARoption.diobs == 0
    % Concatinate all dummies
    Yd = [];
    Xd = [];
elseif BVARoption.socpri == 1 && BVARoption.diobs == 0
    % Concatinate all dummies
    Yd = [Yd1];
    Xd = [Xd1];
elseif BVARoption.socpri == 0 && BVARoption.diobs == 1
    % Concatinate all dummies
    Yd = [Yd2];
    Xd = [Xd2];  
elseif BVARoption.socpri == 1 && BVARoption.diobs == 1
    % Concatinate all dummies
    Yd = [Yd1; Yd2];
    Xd = [Xd1; Xd2];
end

% Number of dummy observations
Td = rows(Yd);

%==========================================================================
% RESULT
%==========================================================================
NDprior.A_prior    = A_prior;
NDprior.a_prior    = a_prior;
NDprior.PHI_prior  = PHI_prior;
NDprior.v_prior    = v_prior;
NDprior.Yd         = Yd;
NDprior.Xd         = Xd;
NDprior.Td         = Td;
