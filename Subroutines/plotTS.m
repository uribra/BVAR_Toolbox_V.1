% Written by:
% Uriel Braham
% uriel.braham@googlemail.com

function plotTS(VARoption,DATA) 

%**************************************************************************
%   FUNCTION TO PLOT TIME SERIES OF ENODENOUS AND EXOGENOUS VARIBALES
%--------------------------------------------------------------------------
%   INPUT: 
%   - Structure 'VARoption'
%   - Structure 'DATA'
%   - Structure 'TSoption'
%--------------------------------------------------------------------------
%   OUTPUT: 
%   - Time series plot of each endogenous and exogenous variable
%--------------------------------------------------------------------------
%
%**************************************************************************

%-----------------------------PRELIMINARY----------------------------------

% Generate Dates
fo_year = VARoption.starty;
lo_year = VARoption.endy;
fo_period = VARoption.startp;
lo_period = VARoption.endp;
frequency = VARoption.frequency;

%
Yendo = DATA.Yendo;
Yexo = DATA.Yexo;
TSoption.nticks         = ceil((VARoption.endy - VARoption.starty)/4) + 1;
nticks = TSoption.nticks;

out =  DatesOption(fo_year,lo_year,frequency,fo_period,lo_period);

nobs = out.nobs;
%[dates, dates_short] = DatesCreate(fo_year,nobs,frequency,fo_period);

if frequency=='q'
    if lo_period==1
    fo = fo_year + 0.00;
    elseif lo_period==2
    fo = fo_year + 0.25;
    elseif lo_period==3
    fo = fo_year + 0.5;
    elseif lo_period==4
    fo = fo_year + 0.75;
    end
end

if frequency=='m'
    if lo_period==1
    fo = fo_year + 0.00;
    elseif lo_period==2
    fo = fo_year + 1/12;
    elseif lo_period==3
    fo = fo_year + 1/12*2;
    elseif lo_period==4
    fo = fo_year + 1/12*3;
    elseif lo_period==5
    fo = fo_year + 1/12*4;
    elseif lo_period==6
    fo = fo_year + 1/12*5;   
    elseif lo_period==7
    fo = fo_year + 1/12*6;   
    elseif lo_period==8
    fo = fo_year + 1/12*7;
    elseif lo_period==9
    fo = fo_year + 1/12*8;
    elseif lo_period==10
    fo = fo_year + 1/12*9;
    elseif lo_period==11
    fo = fo_year + 1/12*10;
    elseif lo_period==12
    fo = fo_year + 1/12*11;    
    end
end


%*************************************************************************%
%                      PLOT ENDOGENOUS VARIABLES                          %
%*************************************************************************%

row_1 = round(sqrt(length(VARoption.names_endo)));
col_1 = ceil(sqrt(length(VARoption.names_endo)));
plot1 = figure;
set(plot1,'name','Time Series Plot: Endogenous Variables');
set(plot1,'DefaultAxesColorOrder',[0 0 0],'DefaultAxesLineStyleOrder','-','DefaultAxesTitleFontWeight','normal');
for jj=1:1:length(VARoption.names_endo)
    subplot(row_1,col_1,jj)
    plot(DATA.Dates, Yendo(:,jj), 'b', 'LineWidth', 1.0)
    title(VARoption.names_endo{jj}, 'FontSize', 12, 'FontWeight','normal')
    hold on
    box on
    grid off
end 

%filename = ['TSplot_endogenous'];
%FigName = [filename];
%set(gcf, 'Color', 'w');
%export_fig(FigName,'-pdf','-png','-painters')
    
filename = ['TSplot_endogenous'];
FigName = [filename];
path = [cd '/Results'];
saveas(gcf, fullfile(path, FigName), 'pdf')

%*************************************************************************%
%                      PLOT EXOGENOUS VARIABLES                           %
%*************************************************************************%
 
if isempty(Yexo)==0
    row_2 = round(sqrt(length(VARoption.names_exo)));
    col_2 = ceil(sqrt(length(VARoption.names_exo)));
    plot2 = figure;
    set(plot2,'name','Time Series Plot: Exogenous Variables');
    set(plot2,'DefaultAxesColorOrder',[0 0 0],'DefaultAxesLineStyleOrder','-','DefaultAxesTitleFontWeight','normal');
    for jj=1:1:length(VARoption.names_exo)
        subplot(row_2,col_2,jj)
        plot(DATA.Dates, Yexo(:,jj), 'b', 'LineWidth', 1.0)
        title(VARoption.names_exo{jj}, 'FontSize', 12, 'FontWeight','normal')
        hold on
        box on
        grid off
    end
    
    %filename = ['TSplot_endogenous'];
    %FigName = [filename];
    %set(gcf, 'Color', 'w');
    %export_fig(FigName,'-pdf','-png','-painters')
    
    filename = ['TSplot_exogenous'];
    FigName = [filename];
    path = [cd '/Results'];
    saveas(gcf, fullfile(path, FigName), 'pdf')
end




%-------------------------------------------------------------------------
%{

%% PLOT ENDOGENOUS VARIABLES
row_1 = round(sqrt(length(VARoption.names_endo)));
col_1 = ceil(sqrt(length(VARoption.names_endo)));

set(figure,'name','Time Series Plot: Endogenous Variables');
%set(0,'DefaultAxesColorOrder',[0 0 0],'DefaultAxesLineStyleOrder','-')
for jj=1:1:length(VARoption.names_endo)
    subplot(row_1,col_1,jj)
    plot(Yendo(:,jj),'Color', [0 0 1], 'LineWidth', 1.0)
    title(VARoption.names_endo{jj}, 'FontSize', 12, 'FontWeight','normal')
    DatesPlot(fo,nobs,nticks,frequency);
    hold on
    box on;
    grid off;
end 

%% PLOT EXOGENOUS VARIABLES
if isempty(Yexo)==0
    row_2 = round(sqrt(length(VARoption.names_exo)));
    col_2 = ceil(sqrt(length(VARoption.names_exo)));
    set(figure,'name','Time Series Plot: Exogenous Variables');
    
    for jj=1:1:length(VARoption.names_exo)
    subplot(row_2,col_2,jj)
    plot(Yexo(:,jj),'Color', [0 0 1], 'LineWidth', 1.0)
    title(VARoption.names_endo{jj}, 'FontSize', 12, 'FontWeight','normal')
    DatesPlot(fo,nobs,nticks,frequency);
    hold on
    box on;
    grid off;
    end  
    
end 

%}
%-------------------------------------------------------------------------







