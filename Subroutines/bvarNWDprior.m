% Written by:
% Uriel Braham
% uriel.braham@googlemail.com


function [NWDprior] = bvarNWDprior(VARoption,BVARoption,DATA)


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
% DUMMIES
%==========================================================================
% Construct dummy for Litterman prior
Yd1 = il1*[diag(sigma.*deltai); 
      zeros(M*(p-1),M)];  
Jp  = diag([1:p].^lambda3);  
Xd1 = il1*[kron(Jp, diag(sigma)) zeros((M*p),q)];

% Construct dummies for prior on covariance matrix of residual;
Yd2 = [diag(sigma)];
Xd2 = [zeros(M,(M*p)) zeros(M,q)];

% Construct dummy for the constant and exogenous variables
Yd3 = il1*il4*[zeros(q,M)];
Xd3 = il1*il4*[zeros(q,M*p) kron(1,eye(q))];  

% Construct dummy for the sum of the coefficients
Yd4 = il5*diag(muYendo0.*deltai);
Xd4 = [kron(ones(1,p),Yd4) zeros(M,q)];

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
Yd5= il6*muYendo0;  
Xd5=[kron(ones(1,p),Yd5) il6*Xbar];

% Generate dummy obervations
if BVARoption.socpri == 0 && BVARoption.diobs == 0
    % Concatinate all dummies
    Yd = [Yd1; Yd2; Yd3];
    Xd = [Xd1; Xd2; Xd3];
elseif BVARoption.socpri == 1 && BVARoption.diobs == 0
    % Concatinate all dummies
    Yd = [Yd1; Yd2; Yd3; Yd4];
    Xd = [Xd1; Xd2; Xd3; Xd4];
elseif BVARoption.socpri == 0 && BVARoption.diobs == 1
    % Concatinate all dummies
    Yd = [Yd1; Yd2; Yd3; Yd5];
    Xd = [Xd1; Xd2; Xd3; Xd5];  
elseif BVARoption.socpri == 1 && BVARoption.diobs == 1
    % Concatinate all dummies
    Yd = [Yd1; Yd2; Yd3; Yd4; Yd5];
    Xd = [Xd1; Xd2; Xd3; Xd4; Xd5];
end

% Number of dummy observations
Td = rows(Yd);

%==========================================================================
% PRIOR
%==========================================================================
% Prior mean of autorgressive coefficicents
A_prior     = (Xd'*Xd)\(Xd'*Yd);
a_prior     = vec(A_prior);
% Prior mean of variance-covariance matrix
SSEd        = (Yd-Xd*A_prior)'*(Yd-Xd*A_prior);
v_prior     = Td-cols(Yd);
SIGMA_prior = SSEd./(v_prior-size(SSEd,2)-1);

%==========================================================================
% RESULT
%==========================================================================
NWDprior.A_prior    = A_prior;
NWDprior.a_prior    = a_prior;
NWDprior.SSEd       = SSEd;
NWDprior.v_prior    = v_prior;
NWDprior.SIGMA_prior= SIGMA_prior;
NWDprior.Yd         = Yd;
NWDprior.Xd         = Xd;
NWDprior.Td         = Td;


