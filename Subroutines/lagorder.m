% Written by:
% Uriel Braham
% uriel.braham@googlemail.com


function [AIChat,SIChat,HQChat]=lagorder(VARoption,DATA,pmax);

%**************************************************************************
% FUNCTION TO COMPUTE INFORMATION CRITERIA FOR p=0,...,pmax
%--------------------------------------------------------------------------
%   INPUT: 
%   - Yendo: (T x M) matrix of endogenous variables
%   - Yexo: (T x N) matrix of exogenous variables
%   - constant: 1=constant, 0=no constant
%   - trend:  1=linear trend, 0=no linear trend
%   - pmax: maximum lag-order
%--------------------------------------------------------------------------
%   OUTPUT: 
%   - AIChat: lag-order selected by Akaike information criterion 
%   - SIChat: lag-order selected by Schwarz information criterion 
%   - HQChat: lag-order selected vy Hannan-Quinn information criterion
%--------------------------------------------------------------------------
%
%**************************************************************************

%----------------------------DATA PREPERATION------------------------------
Yendo = DATA.Yendo;
Yexo    = DATA.Yexo;

% Get initial dimensions
constant    = VARoption.constant;
trend       = VARoption.trend;

[Traw, M] = size(Yendo);
q = constant + trend + size(Yexo,2);

Y1 = Yendo;
Y2 = Yendo;

% Generate lagged Y matrix. This will be part of the X matrix
Ylag=zeros(Traw,M*pmax);
for ii=1:pmax
    Ylag(pmax+1:Traw,(M*(ii-1)+1):M*ii)=Y2(pmax+1-ii:Traw-ii,1:M);
end

% Generate the regressor matrix
if constant ==1
    if isempty(Yexo)==1
        if trend==1
                X1 = [Ylag(pmax+1:Traw,:) ones(Traw-pmax,1) (pmax+1:Traw)' ];
        elseif trend==0 
                X1 = [Ylag(pmax+1:Traw,:) ones(Traw-pmax,1)];  
         end 
    elseif isempty(Yexo)==0
         if trend ==1
                X1 = [Ylag(pmax+1:Traw,:) ones(Traw-pmax,1) ((pmax+1):Traw)' Yexo(pmax+1:Traw,:)];
         elseif trend==0 
                X1 = [Ylag(pmax+1:Traw,:) ones(Traw-pmax,1) Yexo(pmax+1:Traw,:) ]; 
         end
    end 
elseif constant==0    
    if isempty(Yexo)==1
        if trend==1
                X1 = [Ylag(pmax+1:Traw,:) (pmax+1:Traw)'];
        elseif trend==0 
                X1 = [Ylag(pmax+1:Traw,:)];  
        end 
    elseif isempty(Yexo)==0
        if trend ==1
                X1 = [Ylag(pmax+1:Traw,:) ((pmax+1):Traw)' Yexo(pmax+1:Traw,:) ];
        elseif trend==0 
                X1 = [Ylag(pmax+1:Traw,:) Yexo(pmax+1:Traw,:) ]; 
        end
    end  
end

% Delete first pmax-rows
Y1 = Y1(pmax+1:Traw,:);
% Get size of final matrix X
[Traw3 K] = size(X1); 
% Matricies for estimation
Y = Y1;

%--------------------------------------------------------------------------

AICcrit = zeros(pmax,1);
SICcrit = zeros(pmax,1);
HQCcrit = zeros(pmax,1);
% Evaluate criterion for p=0,...,pmax
for jj=0:pmax

    % Regressor matrix for each p=0,...,pmax
	X=X1(:,1: q + jj*M);
    
    % OLS ESTIMATES (= MLE ESTIMATES)
    A_OLS = (X'*X)\(X'*Y);              
    SSE = (Y - X*A_OLS)'*(Y - X*A_OLS);
    % IC have to be evaluated at the MLE (see Kilian and Lütkepohl 2017,
    % p.56).
    SIGMA_MLE = SSE./(Traw-pmax); 
	kjj=length(vec(A_OLS));   

    logL = -((Traw-pmax)/2)*(M*(1+log(2*pi)) + log(det(SIGMA_MLE)));
    
    AICcrit(jj+1,1)=-2*(logL/(Traw-pmax))+2*(kjj/(Traw-pmax));                      % AIC value
    SICcrit(jj+1,1)=-2*(logL/(Traw-pmax))+kjj*log(Traw-pmax)/(Traw-pmax);         % SIC value
    HQCcrit(jj+1,1)=-2*(logL/(Traw-pmax))+kjj*2*log(log(Traw-pmax))/(Traw-pmax);  % HQC value
	
end;
%--------------------------------------------------------------------------

% Rank models for p = 0,1,2,...,pmax
[critmin,critorder]=min(AICcrit);
AIChat=critorder-1;

[critmin,critorder]=min(HQCcrit);
HQChat=critorder-1;

[critmin,critorder]=min(SICcrit);
SIChat=critorder-1;

