% Written by:
% Uriel Braham
% uriel.braham@googlemail.com

function  [result, VAR] = plotIRF(VARoption,VAR,result,irf_record)

%**************************************************************************
%   FUNCTION TO PLOT IMPLUSE REPONSES FUNCTIONS
%--------------------------------------------------------------------------
%   INPUT: 
%   - Structure 'VARoption'
%   - Structure 'VAR'
%   - Structure 'result
%   - Structure 'irf_record''
%--------------------------------------------------------------------------
%   OUTPUT: 
%   - Plot of Impluse Response Functions for each structural shock
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
if ~exist('irf_record')
    error('You need to provide the record of IRFs in "irf_record".');
end

%**************************************************************************
% RETRIEVE INFORMATION FROM STRUCTURES
%**************************************************************************
ihor            = VARoption.ihor;
bandwidthIRF    = VARoption.bandwidthIRF;
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
% SCALE OF IRFs
%**************************************************************************

irfunits            = cell(1,length(VARoption.names_endo));
VARoption.levelcond = ismember(VARoption.units,'rate');

for rr=1:1:length(VARoption.names_endo)
    if VARoption.levelcond(1,rr) == 1
        irfunits{1,rr} = 'Percentage points';
    elseif  VARoption.levelcond(1,rr) == 0
        irfunits{1,rr} = 'Percent';
        irf_record(rr,:) = cellfun(@(x) x*100,irf_record(rr,:),'un',0);
    end
end

%**************************************************************************
% QUANTILES IRFs
%**************************************************************************
irf_quantiles   = cell(M,M);
lower           = ((100-bandwidthIRF)/2)/100; 
upper           = (100 - (100-bandwidthIRF)/2)/100;
median          = 0.5; 
quantiles       = [lower, median , upper];

for kk=1:1:M % Loop over rows
    for ll=1:1:M % Loop over colums
      irf_quantiles{kk,ll} = quantile(irf_record{kk,ll},quantiles)';  
    end
end
result.irf_quantiles = irf_quantiles;

%**************************************************************************
% PLOT ALL SHOCKS IN DIFFERENT FIGURES
%**************************************************************************
shade1 = 0.8*ones(1,3);
shade2 = 0.7*ones(1,3);
%row = round(sqrt(M));
%col = ceil(sqrt(M));
%col=ceil(M^0.5);
%row=ceil(M/ncolumns);
col = 2;
row = ceil(length(VARoption.names_endo)/col);
nshock = M;
nvars = M;

xlabelfontsize = 10;
ylabelfontsize = 10;

for lll = 1:1:nshock  % Loop over shocks (columns)
 
    figname = ['Impulse Response Functions of Shock to: ' VARoption.names_endo{1,lll}];
    set(figure,'name',figname);
    set(gca,'LineWidth',1)
    
    for kkk=1:1:nvars  % Loop over responses (rows)
    
        plot1 = subplot(row,col,kkk);
        
        set(plot1,'XTick',0:xstep_irf:VARoption.ihor);
        set(plot1, 'Xlim', [0 ihor]);
        set(plot1,'DefaultAxesColorOrder',[0 0 0],'DefaultAxesLineStyleOrder','-','DefaultAxesTitleFontWeight','normal')
        minband=min(irf_quantiles{kkk,lll}(:,1));
        maxband=max(irf_quantiles{kkk,lll}(:,3));
        space=maxband-minband;
        Ymin=minband-0.2*space;
        Ymax=maxband+0.2*space;
        set(plot1,'YLim',[Ymin Ymax])
        
        fpatt = [0:ihor ihor:-1:0]';
        fpatt(:,2) = [irf_quantiles{kkk,lll}(:,3); flipud(irf_quantiles{kkk,lll}(:,1))];
        patch(fpatt(:,1),fpatt(:,2),shade1,'EdgeColor',shade2);
        hold on
        plot(0:ihor,irf_quantiles{kkk,lll}(:,2),'k.-');
        %plot(0:ihor,irf_quantiles{kkk,lll}(:,2),'k-');
        hold on
        plot(0:ihor,zeros(ihor+1),'Color', [0 0 0])
        %set(plot1,'fontsize', 10)
        title(names_endo{kkk}, 'FontSize', 11, 'FontWeight','normal')
        xlabel(xlab_irf, 'fontsize', xlabelfontsize)
        ylabel(irfunits{kkk}, 'fontsize', ylabelfontsize)
       
        grid off
        box on 
        
    
    end % End loop over responses (rows)
    
    ax=axes('Units','Normal','Position',[.11 .075 .85 .88],'Visible','off');
    set(get(ax,'Title'),'Visible','on')
    title(['Schock to: ' names_endo{lll}],'FontSize',12,...%'FontName','Times New Roman',
       'FontWeight','normal');
   
    %filename = ['IRF_'];
    %FigName = [filename num2str(lll)];
	%set(gcf, 'Color', 'w');
    %export_fig(FigName,'-pdf','-png','-painters')
    
    filename = ['IRF_'];
    FigName = [filename num2str(lll)];
    path = [cd '/Results'];
    saveas(gcf, fullfile(path, FigName), 'pdf')
   
   
end % End loop over shocks (columns) 


%**************************************************************************
% PLOT ALL SHOCKS IN ONE FIGURE
%**************************************************************************
nshock = M;
nvars = M;
col = nshock;
row = nvars;
spn = 0;
xlabelfontsize = 10;
ylabelfontsize = 10;

figname = ['Impulse Response Functions'];
set(figure,'name',figname);
set(gca,'LineWidth',1)

for kkk = 1:1:nshock % Start loop over shocks (columns)
    
    for lll=1:1:nvars  % Start loop over responses (rows)
        
        spn = spn + 1;
        plot2 = subplot(col,row,spn);
        
        minband=min(irf_quantiles{kkk,lll}(:,1));
        maxband=max(irf_quantiles{kkk,lll}(:,3));
        space=maxband-minband;
        Ymin=minband-0.2*space;
        Ymax=maxband+0.2*space;
        
        set(plot2,'XTick',0:xstep_irf:VARoption.ihor);
        set(plot2,'Xlim',[0 ihor]);
        set(plot2,'YLim',[Ymin Ymax])
        set(plot2,'DefaultAxesColorOrder',[0 0 0],'DefaultAxesLineStyleOrder','-','DefaultAxesTitleFontWeight','normal')
        
        fpatt = [0:ihor ihor:-1:0]';
        fpatt(:,2) = [irf_quantiles{kkk,lll}(:,3); flipud(irf_quantiles{kkk,lll}(:,1))];
        patch(fpatt(:,1),fpatt(:,2),shade1,'EdgeColor',shade2);
        hold on
        plot(0:ihor,irf_quantiles{kkk,lll}(:,2),'k.-');
        %plot(0:ihor,irf_quantiles{kkk,lll}(:,2),'k-');
        hold on
        plot(0:ihor,zeros(ihor+1),'Color', [0 0 0])
        %set(plot2,'fontsize', 10)
        title(names_endo{kkk}, 'FontSize', 11, 'FontWeight','normal')
        xlabel(xlab_irf, 'fontsize', xlabelfontsize)
        ylabel(irfunits{kkk}, 'fontsize', ylabelfontsize)
        
        grid off
        box on 
        
    end % End loop over responses (rows)

end % End loop over shocks (columns)   

ax=axes('Units','Normal','Position',[.11 .075 .85 .88],'Visible','off');
set(get(ax,'Title'),'Visible','on')
title('Shocks ','FontSize',12,...%'FontName','Times New Roman',
    'FontWeight','normal');
% side supertitle
ylabel('Responses ','FontSize',12,.... %'FontName','Times New Roman',
    'FontWeight','normal');
set(get(ax,'Ylabel'),'Visible','on')

%filename = ['IRF_full'];
%FigName = [filename];
%set(gcf, 'Color', 'w');
%export_fig(FigName,'-pdf','-png','-painters')

filename = ['IRF_full'];
FigName = [filename];
path = [cd '/Results'];
saveas(gcf, fullfile(path, FigName), 'pdf');
%--------------------------------------------------------------------------