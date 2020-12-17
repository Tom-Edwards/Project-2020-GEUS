%% Site J perturbation comparisons


DP_1989 = readtable('DensityProfile_KAN-U_1989_metres.csv');
DP89_depths = table2array(DP_1989(:,1));
DP89_densities = table2array(DP_1989(:,2));

% compare with perturbed 

DP_1989_P6 = readtable('Homo6.csv');
DP89_depths_P6 = table2array(DP_1989_P6(:,1));
DP89_densities_P6 = table2array(DP_1989_P6(:,2));

DP_1989_P12 = readtable('12.csv');
DP89_depths_P12 = table2array(DP_1989_P12(:,1));
DP89_densities_P12 = table2array(DP_1989_P12(:,2));

DP_1989_P14 = readtable('14.csv');
DP89_depths_P14 = table2array(DP_1989_P14(:,1));
DP89_densities_P14 = table2array(DP_1989_P14(:,2));

DP_1989_P16 = readtable('16.csv');
DP89_depths_P16 = table2array(DP_1989_P16(:,1));
DP89_densities_P16 = table2array(DP_1989_P16(:,2));

DP_1989_P18 = readtable('18.csv');
DP89_depths_P18 = table2array(DP_1989_P18(:,1));
DP89_densities_P18 = table2array(DP_1989_P18(:,2));

DP_1989_P20 = readtable('20.csv');
DP89_depths_P20 = table2array(DP_1989_P20(:,1));
DP89_densities_P20 = table2array(DP_1989_P20(:,2));

%%

% DP_1989_single = readtable('single.csv');
% DP89_depths_single = table2array(DP_1989_single(:,1));
% DP89_densities_single = table2array(DP_1989_single(:,2));
% 
% plot(DP89_densities,DP89_depths)
% set(gca,'ydir','reverse')
% hold on
% legend('1','2','3','4','5','6')
% ylim([0 40])
% 
% %plot(DP89_densities_P10,DP89_depths_P10)
% 
% plot(DP89_densities_P14,DP89_depths_P14)
% plot(DP89_densities_P16,DP89_depths_P16)
% plot(DP89_densities_P18,DP89_depths_P18)
%plot(DP89_densities_P6,DP89_depths_P6,'--')
%plot(DP89_densities_single,DP89_depths_single)



figure('DefaultLegendFontSize',18);
D1 = plot(DP89_densities,DP89_depths);
D1.LineWidth = 2;

grid on

ax = gca;
box(ax,'on')
ax.LineWidth = 1.5;
ax.FontSize = 18;
set(ax,'Ydir','reverse')
titled = append('Site J with perturbation');
title(titled)
xlabel('Density (kg m^{⁻3})','Interpreter','tex')
ylabel('Depth (m)')
ylim([0 10])

hold on 
% D2 = plot(DP89_densities_P6,DP89_depths_P6,'--');
% D2.LineWidth = 2;
% 
hold on 
D3 = plot(DP89_densities_P20,DP89_depths_P20,'--');
D3.LineWidth = 2;



legend('Site J', ' 20cm thick perturbation', 'other')
set(0,'defaultfigurepapersize',[29.7 16]);
set(0,'defaultfigurepaperposition',[.25 .25 [29.7 16]-0.5]);
set(0,'DefaultTextInterpreter','none');
set(0, 'DefaultFigureUnits', 'centimeters');
set(0, 'DefaultFigurePosition', [.25 .25 [29.7 16]-0.5]);
filename = ('Site J and thick perturbed density profile ');
fpath = './toms/PerbsWD89_13';
print(gcf,fullfile(fpath,filename),'-djpeg');

hold off

%% Homogenous perturbs


DP_1989 = readtable('DensityProfile_KAN-U_1989_metres.csv');
DP89_depths = table2array(DP_1989(:,1));
DP89_densities = table2array(DP_1989(:,2));

% compare with perturbed 

DP_1989_H10 = readtable('Homo10.csv');
DP89_depths_H10 = table2array(DP_1989_H10(:,1));
DP89_densities_H10 = table2array(DP_1989_H10(:,2));

DP_1989_H12 = readtable('Homo12.csv');
DP89_depths_H12 = table2array(DP_1989_H12(:,1));
DP89_densities_H12 = table2array(DP_1989_H12(:,2));

DP_1989_H14 = readtable('Homo14.csv');
DP89_depths_H14 = table2array(DP_1989_H14(:,1));
DP89_densities_H14 = table2array(DP_1989_H14(:,2));

DP_1989_H16 = readtable('Homo16.csv');
DP89_depths_H16 = table2array(DP_1989_H16(:,1));
DP89_densities_H16 = table2array(DP_1989_H16(:,2));

DP_1989_H18 = readtable('Homo18.csv');
DP89_depths_H18 = table2array(DP_1989_H18(:,1));
DP89_densities_H18 = table2array(DP_1989_H18(:,2));

DP_1989_H20 = readtable('Homo20.csv');
DP89_depths_H20 = table2array(DP_1989_H20(:,1));
DP89_densities_H20 = table2array(DP_1989_H20(:,2));

DP_1989_H25 = readtable('Homo25.csv');
DP89_depths_H25 = table2array(DP_1989_H25(:,1));
DP89_densities_H25 = table2array(DP_1989_H25(:,2));

DP_1989_H30 = readtable('Homo30.csv');
DP89_depths_H30 = table2array(DP_1989_P30(:,1));
DP89_densities_H30 = table2array(DP_1989_P30(:,2));

DP_1989_single = readtable('single.csv');
DP89_depths_single = table2array(DP_1989_single(:,1));
DP89_densities_single = table2array(DP_1989_single(:,2));

plot(DP89_densities,DP89_depths)
set(gca,'ydir','reverse')
hold on
%legend('1','2','3','4','5','6')
%plot(DP89_densities_single,DP89_depths_single)
plot(DP89_densities_H10,DP89_depths_H10)
plot(DP89_densities_H12,DP89_depths_H12)
plot(DP89_densities_H14,DP89_depths_H14)
plot(DP89_densities_H16,DP89_depths_H16)
plot(DP89_densities_H18,DP89_depths_H18)
plot(DP89_densities_H20,DP89_depths_H20)
plot(DP89_densities_H25,DP89_depths_H25)
plot(DP89_densities_H30,DP89_depths_H30)
