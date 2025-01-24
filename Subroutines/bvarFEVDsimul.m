% Written by:
% Uriel Braham
% uriel.braham@googlemail.com

function [fevd_record, VAR] = bvarFEVDsimul(VARoption,VAR,result);

%**************************************************************************
%   FUNCTION TO SIMULATE FORECAST ERROR VARIANCE DECOMPOSITIONS OF 
%   IDENTIFIED STRUCTURAL SHOCKS WITH CHOLESKY DECOMOPOSITION
%--------------------------------------------------------------------------
%   INPUT: 
%   - Structure 'VARoption'
%   - Structure 'VAR'
%   - Structure 'result'
%--------------------------------------------------------------------------
%   OUTPUT: 
%   - Structure 'fevd_record'
%   - Structure 'VAR'
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
K           = VAR.K;
q           = VAR.q;
p           = VARoption.p;
ihorfevd    = VARoption.ihorfevd;
nsave       = VARoption.nsave;
A_draws     = result.A_draws;
a_draws     = result.a_draws;
SIGMA_draws = result.SIGMA_draws;

fevd_record = cell(VAR.M,VAR.M);

disp('')
tic;
disp('Posterior simulation for Forecast Error Variance Decomposition')
for ii=1:1:VARoption.nsave
    
    if mod(ii,1000) == 0 % print iterations
        disp([num2str(ii) ' / ' num2str(nsave) ' draws']);
        %toc;
    end
    
    A_draw = squeeze(A_draws(ii,:,:));
    %A_draw = reshape(a_draws(ii,:),VAR.K,VAR.M);
    SIGMA_draw = squeeze(SIGMA_draws(ii,:,:));  
      
    % COMPATION MATRIX FF for Y(t) = FF*Y(t-1) + u(t)
    JJ = zeros(M*p,M);
    JJ(1:M,1:M) = eye(M);
    FF = zeros(K-q,K-q); 
    %FF(1:M,:) = A((q+1):end,:)';
    FF(1:M,:) = A_draw(1:(end-q),:)';
    FF((M+1):end,1:(K-q-M))= eye(M*(p-1));
    % Shock size
    invD = chol(SIGMA_draw)';
    shock = invD;
    % Calculate impluse responses
    IRF_mat = zeros(M,M,ihorfevd);
    IRF_mat(:,:,1) = shock;
    FFjj = eye(M*p);
    for jj=0:ihorfevd+1
        IRF_mat(:,:,jj+1) = JJ'*FFjj*JJ*shock;
        FFjj = FFjj*FF;  
    end 
    
    V = zeros(M,ihorfevd+1);
    V_ALL = nan(M, ihorfevd+1);
    VarDecomp = cell(M,1); 
    % Calucalte the FEVD
    for kk=1:1:M
        for h = 1:ihorfevd+1
            V(:,h)  = diag(squeeze(IRF_mat(:,kk,h))*squeeze(IRF_mat(:,kk,h))');
            V_ALL(:,h) =  diag(squeeze(IRF_mat(:,:,h))*squeeze(IRF_mat(:,:,h))');
        end
            VarDecomp{kk,:} = cumsum(V,2)./cumsum(V_ALL,2)*100;
    end 

    % Distribute the IRFs in the cell structure 'irf_record'
    for kkk=1:1:M % Loop over colums (shocks)
        for lll=1:1:M % Loop over rows (response variables)
            fevd_record{lll,kkk}(ii,:) = VarDecomp{kkk,1}(lll,:);
        end
    end
    
end


