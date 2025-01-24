% Written by:
% Uriel Braham
% uriel.braham@googlemail.com

function [irf_record, VAR] = bvarIRFsimul(VARoption,VAR,result)

%**************************************************************************
%   FUNCTION TO ESTIMATE SIMULATE ORTHOGONALIZED IMPULSE RESPONSE FUNCTIONS 
%   FOR A BAYESIAN VAR USING A CHOLESKI DECOMOPOTION
%--------------------------------------------------------------------------
%   INPUT: 
%   - Structure 'VARoption'
%   - Structure 'VAR'
%--------------------------------------------------------------------------
%   OUTPUT: 
%   Structure 'irf_record' with
%   - IRF.oirf: (MxMxihor+1) matrix with impluse responses 
%   - 
%   - 
%   - 
%--------------------------------------------------------------------------
%
%**************************************************************************

%----------------------------CHECK INPUTS----------------------------------
if ~exist('VARoption')
    error('You need to provide VAR options in "VARoption".');
end
% If there is VARopt get the vnames
names_endo = VARoption.names_endo;
% Check they are not empty
if isempty(names_endo)
    error('Add labels to endogenous variables in "VARoption".');
end
if ~exist('VAR')
    error('You need to provide estmation information in "VAR".');
end
if ~exist('result')
    error('You need to provide estimation results in "results".');
end


%-----------------------------PRELIMINARY----------------------------------
M           = VAR.M;
N           = VAR.N;
K           = VAR.K;
q           = VAR.q;
p           = VARoption.p;
ihor        = VARoption.ihor;
nsave       = VARoption.nsave;
A_draws     = result.A_draws;
a_draws     = result.a_draws;
SIGMA_draws = result.SIGMA_draws;

% IMPLUSE RESPONE FUNCTIONS
irf_record = cell(M,M);

disp('')
tic;
disp('Posterior simulation for (orthogonalized) Impulse Response Functions')

for ii=1:1:VARoption.nsave
    
    if mod(ii,1000) == 0 % print iterations
        disp([num2str(ii) ' / ' num2str(nsave) ' draws']);
        %toc;
    end
    
    A_draw = squeeze(A_draws(ii,:,:));
    %A_draw = reshape(a_draws(ii,:),VAR.K,VAR.M);
    SIGMA_draw = squeeze(SIGMA_draws(ii,:,:));  
    
    % Calculate impluse responses
    [IRF, VAR] = OIRF(VARoption,VAR,A_draw,SIGMA_draw);
    % Distribute the IRFs in the cell structure 'irf_record'
    % Loop over rows
    for kk=1:1:M
        % Loop over colums
        for ll=1:1:M
            irf = squeeze(IRF.oirf(kk,:,1:ihor+1));
            irf_record{kk,ll}(ii,:) = irf(ll,:);
        end
    end
end
