folder = '..\data\encode-net\'; % data file location
nruns = 1; % number of realizations

% constants
COL_N = 1;
COL_NTILDE = 2;
COL_GOUTSTAR = 3;
COL_GOUT = 4;
COL_ASSORT = 5;
COL_ASSORTSTAR = 6;
COL_GSCC = 7;
COL_GSCCSTAR = 8;
COL_GBT = 9;
COL_GBTSTAR = 10;
N_COLS = 10;
FONT_SIZE = 24;
LEGEND_FONT_SIZE = 18;
XMAX = 5.5;

%%
% load file
B = csvread([folder, 'data_dir_11_m2_l0_dmg0.csv'],0,0,[0 0 0 nruns*N_COLS*585-1]);


%%
set(0,'defaulttextinterpreter','latex')


% set up style and generate plots
close all
nametag = '';

mode = 11;
m = 2;
l = 0;
dmg = 0;

name = 'ENCODE ';
npts = 585;


name = [name, nametag];

tag = 1;
prefix = ['data_dir_', num2str(mode), '_m', num2str(m), '_l', num2str(l), '_dmg', num2str(dmg)];
sty = 'r-';

e = B(COL_N:N_COLS:end);
size(e)
e = transpose(reshape(e,npts,nruns));


for i = 2:N_COLS
    
    x = B(i:N_COLS:end);
   
    if i == 2
        x = x - B(COL_GOUTSTAR:N_COLS:end);
    end
    x = transpose(reshape(x,npts,nruns));

        
    figure(i)
    grid on
    set(gcf, 'Position', [100, 100, 840, 660])

    hold on
    plot(mean(e,1), mean(x,1), sty, 'LineWidth', 3)
    hold off
    
    if ismember(i, [COL_GOUTSTAR, COL_GSCCSTAR, COL_GBTSTAR])
        axis([0 XMAX 0 0.35])
        yticks(0:0.05:0.35)
    elseif ismember(i, [COL_ASSORT, COL_ASSORTSTAR])
        axis([0 XMAX -0.3 0.3])
        yticks(-0.3:0.1:0.3)
    elseif i == 11
        axis([0 XMAX -1 1])
        yticks(-1:0.2:1)
    else
        axis([0 XMAX 0 1])
        yticks(0:0.2:1)
    end

    xlabel '$E/N$'
    xticks(0:2:XMAX)

end


figure(11)
grid on
corr = csvread([folder, 'data_dir_', num2str(mode), '_m', num2str(m), '_l0_dmg0_tag0_corrqd.csv']);
set(gcf, 'Position', [100, 100, 840, 660])
hold on
plot(corr(10:end,1), corr(10:end,2), [sty,'-'], 'LineWidth', 3)
plot(corr(10:end,1), corr(10:end,3), sty, 'LineWidth', 3)
hold off

xlabel('$E/N$')
ylabel('$\mathrm{Corr}(d,q)$')
axis([0 XMAX -1 1])


l = legend('ENCODE (in)', 'ENCODE (out)', 'Location', 'Northwest');
l.FontSize = LEGEND_FONT_SIZE;
set(l,'Interpreter','latex');

for i = [2:N_COLS, 11]
    figure(i)
       
    
    if i == COL_NTILDE
        ylabel('$(n_\mathrm{active} - |\mathrm{GOUT^*}|)/N$')
    elseif i == COL_GOUTSTAR
        ylabel('$|\mathrm{GOUT^*}|/N$')
    elseif i == COL_GOUT
        ylabel('$|\mathrm{GOUT}|/N$')
    elseif i == COL_GSCCSTAR
        ylabel('$|\mathrm{GSCC^*}|/N$')
    elseif i == COL_GSCC
        ylabel('$|\mathrm{GSCC}|/N$')
    elseif i == COL_GBTSTAR
        ylabel('$|\mathrm{GBT^*}|/N$')
    elseif i == COL_GBT
        ylabel('$|\mathrm{GBT}|/N$')
    elseif i == COL_ASSORT
        ylabel('in-out assortativity')
    elseif i == COL_ASSORTSTAR
        ylabel('in-out assort. of damaged subnet')
    end
    
    p = 1;
    if ismember(i, [COL_GOUT, COL_GSCC, COL_GBT, COL_NTILDE])
        loc = 'SouthEast';
    elseif ismember(i, [COL_GOUTSTAR, COL_GSCCSTAR, COL_GBTSTAR])
        loc = 'NorthWest';
    else
        loc = 'NorthWest';
    end
    
    if i == 8 || i == 9
        p = 2;
    end
    
    box on
    set(gca, 'FontSize', FONT_SIZE)
    
    if ismember(i, [COL_NTILDE, COL_GOUT, COL_ASSORT])
        l = legend('ENCODE', 'Location', loc);
        l.FontSize = LEGEND_FONT_SIZE;
        set(l,'Interpreter','latex');
    end
        
    set(gca,'TickLabelInterpreter', 'latex');
    set(gca,'XMinorTick','on','YMinorTick','on')

    
    
end

