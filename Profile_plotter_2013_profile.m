
DP_2013 = readtable('DensityProfile_KAN_U_2013_nanmean.csv');

DP13_depths = table2array(DP_2013(:,1));
DP13_densities = table2array(DP_2013(:,2));
DP13_depths_top10m = table2array(DP_2013(1:101,1));
DP13_densities_top10m = table2array(DP_2013(1:101,2));




fpath = './toms/PerbsWD89_13/DPS';

figure('DefaultLegendFontSize',10);
%s1 = stairs(DP12_densities_top10m,DP12_depths_top10m);
s1 = plot(DP13_densities,DP13_depths);

grid on
s1.LineWidth = 2;
grid on
ax = gca;
ax.FontSize = 10;

box(ax,'on')
ax.LineWidth = 2;
ax.FontSize = 16;
set(ax,'Ydir','reverse')
titled = append('Density profile at KanU 2013');
title(titled)
xlabel('Density (kg m^{⁻3})','Interpreter','tex')
ylabel('Depth (m)')
%cm ice slab added to each metre of firn
ylim([0 20])


set(0,'defaultfigurepapersize',[29.7 16]);
set(0,'defaultfigurepaperposition',[.25 .25 [29.7 16]-0.5]);
set(0,'DefaultTextInterpreter','none');
set(0, 'DefaultFigureUnits', 'centimeters');
set(0, 'DefaultFigurePosition', [.25 .25 [29.7 16]-0.5]);
print(gcf,fullfile(fpath,'KanU_2013'),'-djpeg');


%%
