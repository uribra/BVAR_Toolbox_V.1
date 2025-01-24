% Written by:
% Uriel Braham
% Univerity of Bayreuth
% Universitätsstraße 30
% 95445 Bayreuth
% uriel.braham@googlemail.com

function [S, maxeig] = stability(coefmat,M,p,q)

%**************************************************************************
%   FUNCTION TO CHECK FOR STABLITY OF THE LINEAR VAR WITH COMPANION MATRIX
%--------------------------------------------------------------------------
%   INPUT: 
%   - coefmat: (K x M) matrix of coefficients
%   - M: number of endogenous variables
%   - p: lag-order
%   - q: number of deterministic variables 
%--------------------------------------------------------------------------
%   OUTPUT: 
%   - S: S=0 -> VAR system is instable; S=1 -> VAR system is stable.
%--------------------------------------------------------------------------
%
%**************************************************************************

A = coefmat;
K = q + M*p;
FF = zeros(K-q,K-q); 
%FF(1:M,:) = A((q+1):end,:)';
FF(1:M,:) = A(1:end-q,:)';
FF((M+1):end,1:(K-q-M))= eye(M*(p-1));
ee=max(abs(eig(FF)));
S=ee<1;
%if S==1
%    disp("All eigenvalues of the comapanion matrix lie within the unit circle. VAR system is stable.")
%else 
%    disp("At least one eigenvalue of the comapanion matrix lies outside the unit circle. VAR system is instable.")
%end 
maxeig = ee;


