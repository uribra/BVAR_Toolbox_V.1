% Written by:
% Uriel Braham
% Univerity of Bayreuth
% Universit‰tsstraﬂe 30
% 95445 Bayreuth
% uriel.braham@googlemail.com

function IWdraw = IW(scale,v)
%**************************************************************************
%   FUNCTION TO DRAW FROM AN INVERSE WHISHART DISTRIBUTION
%--------------------------------------------------------------------------
%   INPUT: 
%   - scale: (M x M) matrix of sum of squared errors
%   - v: dregrees of freedom
%--------------------------------------------------------------------------
%   OUTPUT: STRUCTURE 'OLS'
%   - IWdraw: (MxM) matrix of draws
%--------------------------------------------------------------------------
%
%**************************************************************************
%scale = nspd(scale);
%inv_scale = inv(scale);
%M = size(inv_scale,1);
%inv_scale_chol = chol(inv_scale);
%Z    = randn(v,M)*inv_scale_chol;
%draw  = inv(Z'*Z);
%IWdraw = draw;

% function [draw]=iwdraw(S,alpha)
% creates a random draw from an inverse Wishart distribution with scale matrix S and degrees of freedom alpha
% inputs:  - matrix 'S': scale matrix for sigma
%          - integer 'alpha': degrees of freedom for sigma
% outputs: - matrix 'draw': random draw from the inverse Wishart distribution

% first obtain a stabilised lower Cholesky factor of S
scale = nspd(scale);
C=chol(scale,'Lower');

% draw the matrix Z of alpha multivariate standard normal vectors
Z=randn(v,size(scale,1));

% obtain the draw
IWdraw=(C/(Z'*Z))*C';
 