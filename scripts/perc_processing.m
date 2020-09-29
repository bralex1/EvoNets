folder = '..\data\main\'; % location of data files

npts = 10000; % number of data points
nruns = 100; % number of realizations

% defining constants
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
XMAX = 14;

%%
% load in data files -- update if using non-default names
B0 = csvread([folder, 'data_dir_0_m2_l0_dmg0.csv']); % random
B9 = csvread([folder, 'data_dir_9_m2_l0_dmg0.csv']); % sens m = 2
B7 = csvread([folder, 'data_dir_7_m2_l0_dmg0.csv']); % hybrid m = 2
B8 = csvread([folder, 'data_dir_8_m2_l0_dmg0.csv']); % struc m = 2


B18 = csvread([folder, 'data_dir_8_m3_l0_dmg0.csv']); % struc m = 3
B19 = csvread([folder, 'data_dir_9_m3_l0_dmg0.csv']); % sens m = 3
B17 = csvread([folder, 'data_dir_7_m3_l0_dmg0.csv']); % hybrid m = 3


%%
set(0,'defaulttextinterpreter','latex')

close all

l = 0;
dmg = 0;
name = '';
% nametag = ['(L = ', num2str(l), ', a = ', num2str(a), ')'];
nametag = '';
label = 1;


% set up labels and styles for plots
for mode = [0,9,8,7]
for m = [2, 3]
    if (mode == 0 && m == 3)
        continue
    end

if mode == 7
    label = 2;
elseif mode == 8
    label = 3;
elseif mode == 9
    label = 4;
end
    


if mode == 0
    name = 'DER ';
elseif mode == 1
    if dmg == 0
        name = ['DCP (M = ', num2str(m), ')'];
    else
        name = ['DCPB (M = ', num2str(m), ')'];
    end
elseif mode == 6
    name = ['DCPB (M = ', num2str(m), ') (annealed)'];
elseif mode == 7
    name = ['AB hybrid (M = ', num2str(m), ')'];
elseif mode == 8
    name = ['AB size (M = ', num2str(m), ')'];
elseif mode == 9
    name = ['AB sensitivity (M = ', num2str(m), ')'];
elseif mode == 5
    if alpha == 0
        name = ['DCP (M = ', num2str(m), ')'];
    else
        name = ['CPDS (M = ', num2str(m), ', \alpha = ', num2str(alpha), ')'];
    end
end   

name = [name, nametag];

tag = 1;
prefix = ['data_dir_', num2str(mode), '_m', num2str(m), '_l', num2str(l), '_dmg', num2str(dmg)];


if mode == 5
    prefix = [prefix, '_alp', num2str(alpha)];
end



sty = '-';
if mode == 0
    sty = 'k-';
    B = B0;
elseif mode == 7
    if m == 3
        sty = 'r--';
        B = B17;
    else
        sty = 'r-';
        B = B7;
    end
elseif mode == 8
    if m == 3
        sty = 'b--';
        B = B18;
    else
        sty = 'b-';
        B = B8;
    end
else
    if m == 3
        sty = 'g--';
        B = B19;
    else
        sty = 'g-';
        B = B9;
    end
end


% format data list as npts x nruns matrix
e = B(COL_N:N_COLS:end);
size(e)
e = transpose(reshape(e,npts,nruns));

% plot figures
for i = 2:N_COLS
    
    x = B(i:N_COLS:end);
   
    if i == 2
        x = x - B(COL_GOUTSTAR:N_COLS:end);
    end
    x = transpose(reshape(x,npts,nruns));

        
    figure(i)
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
corr = csvread([folder, 'data_dir_', num2str(mode), '_m', num2str(m), '_l0_dmg0_tag0_corrqd.csv']);
set(gcf, 'Position', [100, 100, 840, 660])
hold on
plot(corr(10:end,1), corr(10:end,3), sty, 'LineWidth', 3)
hold off

xlabel('$E/N$')
ylabel('$\mathrm{Corr}(d_{\mathrm{out}},q)$')
axis([0 XMAX -1 1])
end
end

% label figures
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
        l = legend('random', 'sensitivity $(m=2)$', 'sensitivity $(m=3)$', 'structure $(m=2)$', ...
            'structure $(m=3)$', 'hybrid $(m=2)$', 'hybrid $(m=3)$', 'Location', loc);
        l.FontSize = LEGEND_FONT_SIZE;
        set(l,'Interpreter','latex');
    end
        
    set(gca,'TickLabelInterpreter', 'latex');
    set(gca,'XMinorTick','on','YMinorTick','on')

    
    
end

