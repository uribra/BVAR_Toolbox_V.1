% Written by:
% Uriel Braham
% uriel.braham@googlemail.com


function [IRF, VAR] = OIRF(VARoption,VAR,A,SIGMA)

%**************************************************************************
%   FUNCTION TO ESTIMATE ORTHOGONALIZED IMPULSE RESPONSE FUNCTIONS USING
%   A CHOLESKI DECOMPOTION
%--------------------------------------------------------------------------
%   INPUT: 
%   - A:
%   - SIGMA:
%   - Structure 'VARoption'
%   - Structure 'VAR'
%--------------------------------------------------------------------------
%   OUTPUT: 
%   Structure 'IRF' with
%   - IRF.oirf: (MxMxihor+1) matrix with impluse responses 
%   Structure 'VAR'
%   - 
%   - 
%   - 
%--------------------------------------------------------------------------
%
%**************************************************************************




% VAR
M           = VAR.M;
N           = VAR.N;
q           = VAR.q; 
K           = VAR.K;

% OPTIONS
p           = VARoption.p;
ihor        = VARoption.ihor;


% COMPATION MATRIX FF for Y(t) = FF*Y(t-1) + u(t)
M = size(A,2);   
JJ = zeros(M*p,M);
JJ(1:M,1:M) = eye(M);
FF = zeros(K-q,K-q); 
%FF(1:M,:) = A((q+1):end,:)';
FF(1:M,:) = A(1:(end-q),:)';
FF((M+1):end,1:(K-q-M))= eye(M*(p-1));

IRF_mat = zeros(M,M,ihor);

invD = chol(SIGMA)';
if VARoption.shocksize == 0
    shock = invD;
elseif VARoption.shocksize == 1
    shock = (invD)*diag(1./diag(invD));
end

% IMPLUSE RESPONSE FUNCTION
IRF_mat(:,:,1) = shock;

FFjj = eye(M*p);
for jj=0:ihor+1
    IRF_mat(:,:,jj+1) = JJ'*FFjj*JJ*shock;
    FFjj = FFjj*FF;  
end 

%----------------------------RESULTS---------------------------------------
IRF.oirf    = IRF_mat;
VAR.comp    = FF;
%--------------------------------------------------------------------------







