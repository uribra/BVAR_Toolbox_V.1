% Written by:
% Uriel Braham
% uriel.braham@googlemail.com

function  [result, VAR] = plotFEVD(VARoption,VAR,result,fevd_record)

%**************************************************************************
%   FUNCTION TO PLOT FORECAST ERROR VARIANCE DECOMOSITIONS
%--------------------------------------------------------------------------
%   INPUT: 
%   - Structure 'VARoption'
%   - Structure 'VAR'
%   - Structure 'result
%   - Structure 'fevd_record_record''
%--------------------------------------------------------------------------
%   OUTPUT: 
%   - Plot of Forecast Error Variance Decomposition for each structural shock
%--------------------------------------------------------------------------
%
%**************************************************************************

%**************************************************************************
% CHECK INPUTS
%**************************************************************************
if ~exist('VARoption')
    error('You need to provide VAR options in "VARoption".');
end
% If there is VARopt get the vnames
names_endo = VARoption.names_endo;
% Check they are not empty
if isempty(names_endo)
    error('Add labels to endogenous variables in "VARoption".');
end
if ~exist('result')
    error('You need to provide estimation results in "results".');
end
if ~exist('VAR')
    error('You need to provide estmation information in "VAR".');
end
if ~exist('fevd_record')
    error('You need to provide the record of FEVDs in "fevd_record".');
end

%**************************************************************************
% RETRIEVE INFORMATION FROM STRUCTURES
%**************************************************************************
bandwidthFEVD   = VARoption.bandwidthFEVD;
ihorfevd        = VARoption.ihorfevd;
M               = VAR.M;
if VARoption.frequency == 'q'
    xlab_irf = 'Quarters';
    xstep_irf = 4;
elseif VARoption.frequency == 'm'
    xlab_irf = 'Months';
    xstep_irf = 12;
elseif VARoption.frequency == 'y'
    xlab_irf = 'Years';
    xstep_irf = 1;
end

%**************************************************************************
% QUANTILES FEVDs
%**************************************************************************
fevd_quantiles  = cell(M,M);
lower           = ((100-bandwidthFEVD)/2)/100; 
upper           = (100 - (100-bandwidthFEVD)/2)/100;
median          = 0.5; 
quantiles       = [lower, median , upper];

for kk=1:1:M % Loop over rows
    for ll=1:1:M % Loop over colums
      fevd_quantiles{kk,ll} = quantile(fevd_record{kk,ll},quantiles)';  
    end
end
result.fevd_quantiles = fevd_quantiles;

%**************************************************************************
% PLOT ALL SHOCKS IN DIFFERENT FIGURES
%**************************************************************************
shade1 = 0.8*ones(1,3);
shade2 = 0.7*ones(1,3);
%row = round(sqrt(M));
%col = ceil(sqrt(M));
%col=ceil(M^0.5);
%row=ceil(M/ncolumns);
xlabelfontsize = 10;
ylabelfontsize = 10;

col = 2;
row = ceil(length(VARoption.names_endo)/col);
nshock = M;
nvars = M;

for lll = 1:1:nshock  % Loop over shocks (columns)
    
    figname = ['Forecast Error Variance Decomposition: ' VARoption.names_endo{1,lll}];
    set(figure,'name',figname);
    set(gca,'LineWidth',1)
    

    for kkk=1:1:nvars      % Loop over responses (rows)
    
        plot1 = subplot(row,col,kkk);
        
        set(plot1,'XTick',0:xstep_irf:ihorfevd);
        set(plot1,'Xlim', [0 ihorfevd]);
        set(plot1,'DefaultAxesColorOrder',[0 0 0],'DefaultAxesLineStyleOrder','-','DefaultAxesTitleFontWeight','normal')
        set(plot1,'YLim',[0 100])
        
        fpatt = [0:ihorfevd ihorfevd:-1:0]';
        fpatt(:,2) = [fevd_quantiles{kkk,lll}(:,3); flipud(fevd_quantiles{kkk,lll}(:,1))];
        patch(fpatt(:,1),fpatt(:,2),shade1,'EdgeColor',shade2);
        hold on
        plot(0:ihorfevd,fevd_quantiles{kkk,lll}(:,2),'k.-');
        %plot(0:ihorfevd,fevd_quantiles{kkk,lll}(:,2),'k-');
        %set(plot1,'fontsize', 10)
        title(names_endo{kkk}, 'FontSize', 11, 'FontWeight','normal')
        xlabel(xlab_irf, 'fontsize', xlabelfontsize)
        ylabel('Percent', 'fontsize', ylabelfontsize)
        
        grid off
        box on 
    
    end % End loop over responses (rows)
    
    ax=axes('Units','Normal','Position',[.11 .075 .85 .88],'Visible','off');
    set(get(ax,'Title'),'Visible','on')
    title(['Schocks to: ' names_endo{lll}],'FontSize',12,...%'FontName','Times New Roman',
       'FontWeight','normal');
   
    %filename = ['FEVD_'];
    %FigName = [filename num2str(lll)];
	%set(gcf, 'Color', 'w');
    %export_fig(FigName,'-pdf','-png','-painters')
    
    filename = ['FEVD_'];
    FigName = [filename num2str(lll)];
    path = [cd '/Results'];
    saveas(gcf, fullfile(path, FigName), 'pdf')
   
end % End loop over shocks (columns) 

%**************************************************************************
% PLOT ALL SHOCKS IN ONE FIGURE
%*************************************************************************
nshock = M;
nvars = M;
col = nshock;
row = nvars;
spn = 0;
xlabelfontsize = 10;
ylabelfontsize = 10;

figname = 'Forecast Error Variance Decomposition';
set(figure,'name',figname);
set(gca,'LineWidth',1)


for kkk = 1:1:nshock % Start loop over shocks (columns)
    
    for lll=1:1:nvars     % Start loop over responses (rows)
        
        spn = spn + 1;
        plot2 = subplot(col,row,spn);
        
        set(plot2,'XTick',0:xstep_irf:ihorfevd);
        set(plot2, 'Xlim', [0 ihorfevd]);
        set(plot2,'DefaultAxesColorOrder',[0 0 0],'DefaultAxesLineStyleOrder','-','DefaultAxesTitleFontWeight','normal')
        set(plot2,'YLim',[0 100])
        
        fpatt = [0:ihorfevd ihorfevd:-1:0]';
        fpatt(:,2) = [fevd_quantiles{kkk,lll}(:,3); flipud(fevd_quantiles{kkk,lll}(:,1))];
        patch(fpatt(:,1),fpatt(:,2),shade1,'EdgeColor',shade2);
        hold on
        plot(0:ihorfevd,fevd_quantiles{kkk,lll}(:,2),'k.-');
        %plot(0:ihor,irf_quantiles{kkk,lll}(:,2),'k-');
        %set(plot2,'fontsize', 10)
        title(names_endo{kkk}, 'FontSize', 11, 'FontWeight','normal')
        xlabel(xlab_irf, 'fontsize', xlabelfontsize)
        ylabel('Percent', 'fontsize', ylabelfontsize)
        
        grid off
        box on 
        
    end    % End loop over responses (rows)
end % End loop over shocks (columns)

ax=axes('Units','Normal','Position',[.11 .075 .85 .88],'Visible','off');
set(get(ax,'Title'),'Visible','on')
title('Shocks ','FontSize',12,...%'FontName','Times New Roman',
    'FontWeight','normal');
% side supertitle
ylabel('Responses ','FontSize',12,.... %'FontName','Times New Roman',
    'FontWeight','normal');
set(get(ax,'Ylabel'),'Visible','on')

%filename = ['FEVD_full'];
%FigName = [filename];
%set(gcf, 'Color', 'w');
%export_fig(FigName,'-pdf','-png','-painters')

filename = ['FEVD_full'];
FigName = [filename];
path = [cd '/Results'];
saveas(gcf, fullfile(path, FigName), 'pdf');
%--------------------------------------------------------------------------