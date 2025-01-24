% Written by:
% Uriel Braham
% uriel.braham@googlemail.com

function plotFIT(VARoption,DATA,result) 

%**************************************************************************
%   FUNCTION TO PLOT FITTED VALUES OF ESTIMATED THE BVAR MODEL
%--------------------------------------------------------------------------
%   INPUT: 
%   - Structure 'VARoption'
%   - Structure 'TSoption'
%   - Structure 'DATA'
%   - Structure 'result'
%--------------------------------------------------------------------------
%   OUTPUT: 
%   - Plot of fitted and actucal values of endogenous varibales in the BAVR
%--------------------------------------------------------------------------
%
%**************************************************************************

%----------------------------CHECK INPUTS----------------------------------
if ~exist('VARoption')
    error('You need to provide VAR options in "VARoption".');
end
if ~exist('DATA')
    error('You need to provide the the data in "DATA".');
end
if ~exist('result')
    error('You need to provide estimation results in "results".');
end

%-----------------------------PRELIMINARY----------------------------------
% Unpack necessary stuff
p       = VARoption.p;
Yendo   = DATA.Yendo;
FITTED  = result.FITTED;

% Generate Dates
fo_year = VARoption.starty;
lo_year = VARoption.endy;
fo_period = VARoption.startp;
lo_period = VARoption.endp;
frequency = VARoption.frequency;

%

nticks         = ceil((VARoption.endy - VARoption.starty)/4) + 1;
%nticks = TSoption.nticks;
%[dates, dates_short] = DatesCreate(fo_year,nobs,frequency,fo_period);


if lo_period==1
    fo = fo_year + 0.00;
elseif lo_period==2
    fo = fo_year + 0.25;
elseif lo_period==3
    fo = fo_year + 0.5;
elseif lo_period==4
    fo = fo_year + 0.75;
end


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

if frequency=='q'
    fo1 = fo + VARoption.p*(1/4);
    fo_year = floor(fo1 - VARoption.p*(1/4));
    fo_period = (fo + VARoption.p*(1/4) - fo_year) - 1;
elseif frequency=='m'
    fo1 = fo + VARoption.p*(1/12) + 1/12;
    fo_year = floor(fo1 - VARoption.p*(1/12));
    fo_period = (fo + VARoption.p*(1/12) - fo_year) - 1;
end

out =  DatesOption(fo_year,lo_year,frequency,fo_period,lo_period);
nobs = out.nobs;
%*************************************************************************%
%                    PLOT FITTED AND ACTUAL VALUES                        %
%*************************************************************************%

col = 2;
row = ceil(length(VARoption.names_endo)/col);
plot1 = figure;
set(plot1,'name','Fitted vs Actual Values');
set(plot1,'DefaultAxesColorOrder',[0 0 0],'DefaultAxesLineStyleOrder','-','DefaultAxesTitleFontWeight','normal');
for jj=1:1:length(VARoption.names_endo)
    subplot(row,col,jj)
    plot(DATA.Dates((p+1):end,1), Yendo((p+1):end,jj), 'b', ...
    DATA.Dates((p+1):end,1), FITTED(:,jj),'r', 'LineWidth', 1.0)
    title(VARoption.names_endo{jj}, 'FontSize', 12, 'FontWeight','normal')
    box on
    grid off
    lgd = legend('Actual','Fitted','Location','Best');
    lgd.Box = 'off';
end 

%filename = ['fittedVactual'];
%FigName = [filename];
%set(gcf, 'Color', 'w');
%export_fig(FigName,'-pdf','-png','-painters')
    
filename = ['fittedVactual'];
FigName = [filename];
path = [cd '/Results'];
saveas(gcf, fullfile(path, FigName), 'pdf')

%-------------------------------------------------------------------------%
%{
row = round(sqrt(length(VARoption.names_endo)));
col = ceil(sqrt(length(VARoption.names_endo)));
set(figure,'name','Fitted vs Actual Values');
for jj=1:1:length(VARoption.names_endo)
    subplot(row,col,jj)
    plot([Yendo((VARoption.p+1):end,jj)],'b', 'LineWidth', 1.0)
    hold on
    plot([FITTED(:,jj)],'r', 'LineWidth', 1.0)
    title(VARoption.names_endo{jj}, 'FontSize', 12, 'FontWeight','normal')
    DatesPlot(fo1,nobs,nticks,frequency);
    box on
    grid off
    lgd = legend('Actual','Fitted','Location','Best');
    lgd.Box = 'off';
end
%}
%-------------------------------------------------------------------------%







