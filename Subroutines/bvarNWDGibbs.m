% Written by:
% Uriel Braham
% uriel.braham@googlemail.com

function [result, VAR] = bvarNWDGibbs(VARoption,BVARoption,VAR,result)

%**************************************************************************
%   FUNCTION TO DRAW FROM POSTERIOR DISTIBUTIONS WITH GIBBS SAMPLER
%--------------------------------------------------------------------------
%   INPUT: Structure 'DATA' with
%   -
%--------------------------------------------------------------------------
%   OUTPUT: STRUCTURE OF NWD RESULTS
%   -   
%--------------------------------------------------------------------------
%
%**************************************************************************

%*************************************************************************%
%                UNPACK structures 'VARoption' and 'VAR'                  %
%*************************************************************************%
constant        = VARoption.constant;
p               = VARoption.p;
trend           = VARoption.trend;
nsave           = VARoption.nsave;
nburn           = VARoption.nburn;
ntot            = nburn + nsave;
M               = VAR.M;
K               = VAR.K;
k               = VAR.k;
q               = VAR.q;
N               = VAR.N;
T               = VAR.T; 

%*************************************************************************%
%                             POSTERIOR                                   %
%*************************************************************************%
A_post          = result.A_post;
a_post          = result.a_post; 
SCALE_post      = result.SCALE_post;
v_post          = result.v_post;
SIGMA_post      = result.SIGMA_post; 

Xs              = VAR.Xs;
Ys              = VAR.Ys;
XsXs            = Xs'*Xs;
invXsXs         = XsXs\eye(cols(XsXs));
Ts              = VAR.Ts;
%v_post          = Ts - K + 2;

%*************************************************************************%
%                             GIBBS SAMPLER                               %
%*************************************************************************%
% Storage space for posterior draws
a_draws         = zeros(nsave,K*M);   
A_draws         = zeros(nsave,K,M);  
SIGMA_draws     = zeros(nsave,M,M);

% OLS Estimates
Yreg            = VAR.Yreg;
Xreg            = VAR.Xreg;
A_OLS           = (Xreg'*Xreg)\(Xreg'*Yreg);
SIGMA_OLS       = ((Yreg-Xreg*A_OLS)'*(Yreg-Xreg*A_OLS))./(T-K);
% Initial value of the Gibbs sampler
SIGMA_draw      = SIGMA_OLS; 


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
    % 1.) Draw a|SIGMA,y ~ N(a_post, kron(SIGMA,invXXs))
    V_post = kron(SIGMA_draw,invXsXs);
    a_draw  = a_post + (randn(1,M*K)*chol(V_post))';
    A_draw = reshape(a_draw,K,M);
    %U  = Ys - Xs*A_draw;
    
    % 2.) Draw SIGMA|y ~ iW(scale,v_post)
    %scale  = U'*U;
    SIGMA_draw = IW(SCALE_post,v_post);
    %======================================================================
    
    if irep > nburn        
    %Save draws of the parameters
    a_draws(irep-nburn,:)       = a_draw;
    A_draws(irep-nburn,:,:)     = A_draw;
    SIGMA_draws(irep-nburn,:,:) = SIGMA_draw; 
    % End saving  
    end
% End Gibbs loop
end

A_SD     = squeeze(std(A_draws,1)); 
SIGMA_SD = squeeze(std(SIGMA_draws,1));

%*************************************************************************%
%                             RESULTS                                     %
%*************************************************************************%
result.A_draws     = A_draws;
result.a_draws     = a_draws;
result.SIGMA_draws = SIGMA_draws;
result.A_postSD    = A_SD;
result.SIGMA_postSD=SIGMA_SD;


