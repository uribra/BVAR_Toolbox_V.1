% Written by:
% Uriel Braham
% uriel.braham@googlemail.com

function [VAR, DATA] = VARgen(VARoption,DATA)

constant    = VARoption.constant;              
trend       = VARoption.trend;                 
p           = VARoption.p;   

% Endogenous variables
position_index_endo = nan(1, size(VARoption.names_endo, 2));
for jj=1:1:size(VARoption.names_endo, 2)
    variable_name_endo = VARoption.names_endo(1,jj);
    position_index_endo(1,jj)= find(strcmp(DATA.Series,variable_name_endo));  
end

Yendo       = DATA.VARS(:,[position_index_endo]);
M           = length(VARoption.names_endo);



% Exoegenous variables
if isempty(VARoption.names_exo)==0
    position_index_exo = nan(1, size(VARoption.names_exo, 2));
    for jj=1:1:size(VARoption.names_endo, 2)
        variable_name_exo = VARoption.names_exo(1,jj);
        position_index_exo(1,jj)= find(strcmp(DATA.Series,variable_name_exo));  
    end
    Yexo     = DATA.VARS(:,position_index_exo);
else 
    Yexo     = [];
end
N           = length(VARoption.names_exo);

Traw        = rows(Yendo); 
T           = Traw -p;

% Total number of exogenous variables 
q           = constant + trend + N;
% Coefficents per equation 
K           = M*p + q;
% Total number of coefficients
k           = M*K;

VAR.M       = M;
VAR.N       = N;
VAR.q       = q;
VAR.k       = k;
VAR.T       = T;
DATA.Yendo  = Yendo;
DATA.Yexo   = Yexo;








               