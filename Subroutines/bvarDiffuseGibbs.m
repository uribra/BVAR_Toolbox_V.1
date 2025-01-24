% Written by:
% Uriel Braham
% uriel.braham@googlemail.com

function [result, VAR] = bvarDiffuseGibbs(VARoption,VAR,result)

%**************************************************************************
%   FUNCTION TO SIMULATE (DRAW) FROM THE POSTERIOR DISTRIBUTIONS FOR THE
%   COEFFICIENTS
%--------------------------------------------------------------------------
%   INPUT: 
%   - Structure 'VARoption'
%   - Structure 'VAR'
%   - STructure 'result'
%--------------------------------------------------------------------------
%   OUTPUT: 
%   - Structure 'result' with estimation result of the BAVR
%   - Structure 'VAR' with estimation information
%--------------------------------------------------------------------------
%
%**************************************************************************

%----------------------------CHECK INPUTS----------------------------------
if ~exist('VARoption')
    error('You need to provide BVAR options in "VARoption".');
end
if ~exist('VAR')
    error('You need to provide estmation information in "VAR".');
end
if ~exist('result')
    error('You need to provide estimation results in "results".');
end
%--------------------------------------------------------------------------

%*************************************************************************%
%                UNPACK structures 'VARoption' and 'VAR'                  %
%*************************************************************************%
% Gibbs sampler options
nsave       = VARoption.nsave;
nburn       = VARoption.nburn;
ntot        = nburn + nsave;
it_print    = 2000;
% VAR information
M           = VAR.M;
K           = VAR.K;
k           = VAR.k;
T           = VAR.T;
p           = VARoption.p;
X           = VAR.Xreg;
Y           = VAR.Yreg;
% Posteriors
A_post      = result.A_post;
a_post      = result.a_post;
SSE         = result.SSE;    
SIGMA_post  = result.SIGMA_post; 
v_post      = T-K;
% Storage space for posterior draws
a_draws     = zeros(nsave,K*M);   
A_draws     = zeros(nsave,K,M);  
SIGMA_draws = zeros(nsave,M,M);  

%*************************************************************************%
%                             GIBBS SAMPLER                               %
%*************************************************************************%
% Initial value: 
SIGMA_draw = SIGMA_post;

disp('')
disp('Posterior simulation for BVAR coefficients (Gibbs sampler):')
tic;
%Start Gibbs "loop"
for irep = 1:ntot  % loop over irep
    
    if mod(irep,1000) == 0 % print iterations
        disp([num2str(irep) ' / ' num2str(ntot) ' draws']);
        %toc;
    end
    
    %======================== POSTERIOR DRAWS =============================
    % 1.) Posterior of a|SIGMA,Y ~ N(a_post,V_post) with a_post = a_OLS and
    % V_post = kron(SIGMA_draw, (X'*X)^(-1))
    % Draw a
    V_post = kron(SIGMA_draw,inv(X'*X));
    a_draw = a_post + chol(V_post)'*randn(K*M,1);
    A_draw = reshape(a_draw,K,M);
        
    % 2.) Posterior of SIGMA|Y ~ iW(SSE,T-K) with SSE =(Y-X*A_post)'*(Y-X*A_post)
    % and v_post = T-K;
    % Draw SIGMA
    SIGMA_draw = IW(SSE,v_post);
    %======================================================================
    
    % Store results  
    if irep > nburn  
        a_draws(irep-nburn,:)       = a_draw;
        A_draws(irep-nburn,:,:)     = A_draw;
        SIGMA_draws(irep-nburn,:,:) = SIGMA_draw; 
    end % End Storage    
end % End Gibbs "loop"    

A_SD     = squeeze(std(A_draws,1)); 
SIGMA_SD = squeeze(std(SIGMA_draws,1));

%*************************************************************************%
%                             RESULTS                                     %
%*************************************************************************%
result.A_draws     = A_draws;
result.a_draws     = a_draws;
result.A_postSD    = A_SD;
result.SIGMA_draws = SIGMA_draws;
result.SIGMA_postSD=SIGMA_SD;
