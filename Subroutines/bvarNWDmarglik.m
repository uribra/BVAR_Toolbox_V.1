% Written by:
% Uriel Braham
% uriel.braham@googlemail.com


function [lnmarglik] = bvarNWDmarglik(Xd,Yd,Xs,Ys)

Td = rows(Yd);
Ts = rows(Ys);
M = cols(Yd);
K = cols(Xd);

Phi0   = (inv(Xd'*Xd))*(Xd'*Yd);
S0     = (Yd'*Yd)-Yd'*Xd*(inv(Xd'*Xd))*Xd'*Yd;

Phi1   = (inv(Xs'*Xs))*(Xs'*Ys); 
S1     = (Ys'*Ys)-Ys'*Xs*(inv(Xs'*Xs))*Xs'*Ys;

%* compute constants for integrals

i=1;  
gam0=0; 
gam1=0; 

while i <= M
   gam0=gam0+log(gamma(0.5*(Td-K+1-i))); 
   gam1=gam1+log(gamma(0.5*(Ts-K+1-i))); 
   i=i+1; 
end

%** dummy observation
lnpY0 = -M*(Td-K)*0.5*log(pi)-(M/2)*log(abs(det(Xd'*Xd)))-(Td-K)*0.5*log(abs(det(S0)))+M*(M-1)*0.25*log(pi)+gam0;

%** dummy and actual observations   
lnpY1 = -M*(Ts-K)*0.5*log(pi)-(M/2)*log(abs(det(Xs'*Xs)))-(Ts-K)*0.5*log(abs(det(S1)))+M*(M-1)*0.25*log(pi)+gam1;  
lnmarglik = lnpY1-lnpY0;





