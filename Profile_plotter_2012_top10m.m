
DP_2012 = readtable('RetMIP_density_KAN-U 2012.csv');

DP12_depths = table2array(DP_2012(:,1));
DP12_densities = table2array(DP_2012(:,2));
DP12_depths_top10m = table2array(DP_2012(1:101,1));
DP12_densities_top10m = table2array(DP_2012(1:101,2));

fpath = './toms/PerbsWD89_13/DPS';

figure('DefaultLegendFontSize',10);
s1 =  stairs(DP12_densities_top10m,DP12_depths_top10m);
%s1 = stairs(DP12_densities,DP12_depths);
s1.LineWidth = 2;
grid on
ax = gca;
ax.FontSize = 10;

box(ax,'on')
ax.LineWidth = 2;
ax.FontSize = 16;
set(ax,'Ydir','reverse')
titled = append('Composite Density profile KanU 2012 & SiteJ 1989');
title(titled)
xlabel('Density (kg m^{⁻3})','Interpreter','tex')
ylabel('Depth (m)')
%cm ice slab added to each metre of firn
ylim([0 10.1])


set(0,'defaultfigurepapersize',[29.7 16]);
set(0,'defaultfigurepaperposition',[.25 .25 [29.7 16]-0.5]);
set(0,'DefaultTextInterpreter','none');
set(0, 'DefaultFigureUnits', 'centimeters');
set(0, 'DefaultFigurePosition', [.25 .25 [29.7 16]-0.5]);
print(gcf,fullfile(fpath,'KanU_2012_comp top 10m'),'-djpeg');


%%
