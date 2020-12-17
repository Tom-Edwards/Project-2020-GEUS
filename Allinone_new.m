%% NetCDF file reader
% Homo goes from 6:2:20 + 25, 30,35
% Thick goes through 10:2:30
% Remember to change the titles of the graphs when changing perturbation
% type

number = PT1;      % Perturbation number
Type = PTT;        % Perturbation type
number = num2str(number);
pert = Type;       % Thick or Homo, (slab addition or homogeneous perturbation)
heads = append(pert,number);
outputFP = append('/PerturbationsWD2012_1989/P',heads);
ncfile_firn_density = append('./Output',outputFP,'/rho_bin_1.nc');
ncfile_ice_content = append('./Output',outputFP,'/snic_bin_1.nc');
ncfile_snow_content = append('./Output',outputFP,'/snowc_bin_1.nc');
rho = ncread(ncfile_firn_density,'rho');
Depth = ncread(ncfile_firn_density,'Depth');
ice = ncread(ncfile_ice_content,'snic');
snow = ncread(ncfile_snow_content,'snowc');
DP_2012 = readtable('RetMIP_density_KAN-U 2012.csv');
DP_2013 = readtable('DensityProfile_KAN_U_2013_nanmean.csv');
DP_1989 = readtable('DensityProfile_KAN-U_1989_metres.csv');

LW = 2; % Line width for all plots
Year = '2012';  % Changeable
SY13 = 212575;  % 2012 1st May 204525
SY12 = 204525;  % 2013 + 8050 = 212575
SY = SY12;      % Sets year of final simulation

fpath = append('./toms/PerbsWD89_13/P',heads);
Txtsz = 20;

%Colors

c = [118/255 133/255 212/255];  % blue 
d = [243/255 119/255 56/255];   % orange
e = [0.5020    0.5020    0.5020];  % grey

%% Calculation for bulk density. 

rho_all = (snow + ice)./(snow./rho + ice/917);

%% Variable generator

depth_layers_final = Depth(:,SY);
depth_layers_init = Depth(:,1);
rho_sim_init = rho_all(:,1);
rho_sim_final = rho_all(:,SY);

% 2013 data
x1_13 = DP_2013{:,2};
y1_13 = DP_2013{:,1};
x1_12 = DP_2012{:,2};
y1_12 = DP_2012{:,1};
x1_89 = DP_1989{1:10000,2};
y1_89 = DP_1989{1:10000,1};

%% 2012 DATA

AVGD_per_metre_2012 = table2array(readtable('2012AVGD_per_metre.txt'));
AVGD_per_halfmetre_2012 = table2array(readtable('2012AVGD_per_halfmetre.txt'));
AVGD_per_2metres_2012 = table2array(readtable('2012AVGD_per_2metres.txt'));
AVGD_per_5metres_2012 = table2array(readtable('2012AVGD_per_5metres.txt'));
AVGD_per_10cm_2012 = table2array(readtable('2012AVGD_per_10cm.txt'));

%% 2012 Data regional top 10m

AVGD_2012_Top10m_20cm_RES = table2array(readtable('2012 average density per 20 centimetres.txt'));
AVGD_2012_Top10m_hfm_RES = table2array(readtable('2012 average density per half metre.txt'));
AVGD_2012_Top10m_1m_RES = table2array(readtable('2012 average density per metre.txt'));

%% 2013 Data


AVGD_per_metre_2013 = table2array(readtable('2013 average density per metre.txt'));
AVGD_per_halfmetre_2013 = table2array(readtable('2013 average density per halfmetre.txt'));
AVGD_per_2metres_2013 = table2array(readtable('2013 average density per 2metres.txt'));
AVGD_per_5metres_2013 = table2array(readtable('2013 average density per 5metres.txt'));

AVGD_per_metre_Top5m_2013 = table2array(readtable('RMSD_TOP5m_metre_res.txt'));
AVGD_per_hfmetre_Top5m_2013 = table2array(readtable('RMSD_TOP5m_hf_metre_res.txt'));
AVGD_per_20cm_Top5m_2013 = table2array(readtable('RMSD_TOP5m_twenty_cm_res.txt'));

AVGD_per_metre_5to10m_2013 = table2array(readtable('RMSD_5to10m_metre_res.txt'));
AVGD_per_hfmetre_5to10m_2013 = table2array(readtable('RMSD_5to10m_hf_metre_res.txt'));
AVGD_per_20cm_5to10m_2013 = table2array(readtable('RMSD_5to10m_twenty_cm_res.txt'));

AVGD_per_metre_10to15m_2013 = table2array(readtable('RMSD_10to15m_metre_res.txt'));
AVGD_per_hfmetre_10to15m_2013 = table2array(readtable('RMSD_10to15m_hf_metre_res.txt'));
AVGD_per_20cm_10to15m_2013 = table2array(readtable('RMSD_10to15m_twenty_cm_res.txt'));

AVGD_per_metre_15to19m_2013 = table2array(readtable('RMSD_2013_15m_to_19m_metre_res.txt'));
AVGD_per_hfmetre_15to19m_2013 = table2array(readtable('RMSD_2013_15m_to_19m_hf_metre_res.txt'));
AVGD_per_20cm_15to19m_2013 = table2array(readtable('RMSD_2013_15m_to_19m_twenty_cm_res.txt'));

%% Interpolator final time-step sim values - Take top 20m of firn, round to nearest cm and interpolate

Top20_Depth_final = depth_layers_final(1:90);
Rounded_Depth_final= round(Top20_Depth_final,2);
y = linspace(min(Rounded_Depth_final),max(Rounded_Depth_final),200);
rounded_again_final = round(y,2);
rounded_Depthfinal = rounded_again_final';

Top20_Densities_final = rho_sim_final(1:90);

x_Depths_Final = Rounded_Depth_final;
y_Densities_Final = Top20_Densities_final;
xi_interpolated_depths_final = rounded_Depthfinal;
yi_interpolated_densities_final = interp1(x_Depths_Final,y_Densities_Final,xi_interpolated_depths_final,'next');

%% Interpolator final time-step sim values - Take top 20m of firn, round to nearest cm and interpolate

Top20_Depth_init = depth_layers_init(1:90);
Rounded_Depth_init= round(Top20_Depth_init,2);
y = linspace(min(Rounded_Depth_init),max(Rounded_Depth_init),200);
rounded_again_init = round(y,2);
rounded_Depthinit = rounded_again_init';

Top20_Densities_init = rho_sim_init(1:90);

x_Depths_init = Rounded_Depth_init;
y_Densities_init = Top20_Densities_init;
xi_interpolated_depths_init = rounded_Depthinit;
yi_interpolated_densities_init = interp1(x_Depths_init,y_Densities_init,xi_interpolated_depths_init,'next');


%% Discretizing and RMS
% simulation contains the variables depth_layers_final, ind (resolution),
% rho_sim_final, rho_sim_init, 'res' sets resolution of avg dens per 'r'
 
% indices = discretize(xi_interpolated_depths_final,res);
% good_ind = 1:(find(isnan(indices),1,'first')-1);
% rho_simulated_finalDP = accumarray(indices(good_ind),yi_interpolated_densities_final(good_ind),[],@mean);

%% This does all of what happens above but in less code via the function

% Selects only the top 20m of firn, discretises and interpolates depth and  
% density simultaneously. 

% TO DO - 'r' can be a string in the file name i.e simulated_finalrho_AVGDP'r'M
% sometimes there can be a zero entry as the depth values might skip some
% metres in the deeper parts of the firn. 

%% Ten centimetre resolution

P = 0.1;
res = 0:P:20;

simulated_finalrho_AVGDP10cm = ConvertDepth(xi_interpolated_depths_final,yi_interpolated_densities_final,res);
simulated_initialrho_AVGDP10cm = ConvertDepth(xi_interpolated_depths_init,yi_interpolated_densities_init,res);

%% ten centimetre resolution

P = 0.2;
res = 0:P:20;

simulated_finalrho_AVGDP20cm = ConvertDepth(xi_interpolated_depths_final,yi_interpolated_densities_final,res);
simulated_initialrho_AVGDP20cm = ConvertDepth(xi_interpolated_depths_init,yi_interpolated_densities_init,res);

%% half metre resolution

P = 0.5;
res = 0:P:20;

simulated_finalrho_AVGDPHM = ConvertDepth(xi_interpolated_depths_final,yi_interpolated_densities_final,res);
simulated_initialrho_AVGDPHM = ConvertDepth(xi_interpolated_depths_init,yi_interpolated_densities_init,res);

%% 1 metre resolution

P = 1;
res = 0:P:20;

simulated_finalrho_AVGDP1M = ConvertDepth(xi_interpolated_depths_final,yi_interpolated_densities_final,res);
simulated_initialrho_AVGDP1M = ConvertDepth(xi_interpolated_depths_init,yi_interpolated_densities_init,res);

A1M = numel(simulated_initialrho_AVGDP1M);
B1M = numel(simulated_finalrho_AVGDP1M);
%filepath = {'./python_outputs/Site_J_perturbations_iteration'};

%% 2 metre resolution

P = 2;
res = 0:P:20;

simulated_finalrho_AVGDP2M = ConvertDepth(xi_interpolated_depths_final,yi_interpolated_densities_final,res);
simulated_initialrho_AVGDP2M = ConvertDepth(xi_interpolated_depths_init,yi_interpolated_densities_init,res);

%% 5 metre resolution

P = 5;
res = 0:P:20;

simulated_finalrho_AVGDP5M = ConvertDepth(xi_interpolated_depths_final,yi_interpolated_densities_final,res);
simulated_initialrho_AVGDP5M = ConvertDepth(xi_interpolated_depths_init,yi_interpolated_densities_init,res);

%% RMSD between final and initial simulated profiles

%RMSD_per10cm_sims = mean((simulated_finalrho_AVGDP10cm-simulated_initialrho_AVGDP10cm).^2);
RMSD_perhalfmeter_sims = sqrt(mean((simulated_finalrho_AVGDPHM-simulated_initialrho_AVGDPHM).^2));
RMSD_permeter_sims = sqrt(mean((simulated_finalrho_AVGDP1M-simulated_initialrho_AVGDP1M).^2));
RMSD_per2meters_sims = sqrt(mean((simulated_finalrho_AVGDP2M-simulated_initialrho_AVGDP2M).^2));
RMSD_per5meters_sims = sqrt(mean((simulated_finalrho_AVGDP5M-simulated_initialrho_AVGDP5M).^2));

RMSD_final_initial_sims = table(RMSD_perhalfmeter_sims,...
    RMSD_permeter_sims,RMSD_per2meters_sims,RMSD_per5meters_sims);
filename_append = append(Year,' RMSD_Fin_init_sims_sqrt');
writetable(RMSD_final_initial_sims,fullfile(fpath,filename_append));

%% RMSD between final simulated profile and 2012 profile at KANU

%RMSD_per_10centimeter_2012_finalsim = mean((simulated_finalrho_AVGDP10cm-AVGD_per_10cm_2012).^2);
RMSD_per_halfmeter_2012_finalsim = sqrt(mean((simulated_finalrho_AVGDPHM-AVGD_per_halfmetre_2012).^2));
RMSD_per_meter_2012_finalsim = sqrt(mean((simulated_finalrho_AVGDP1M-AVGD_per_metre_2012).^2));
RMSD_per_2meters_2012_finalsim = sqrt(mean((simulated_finalrho_AVGDP2M-AVGD_per_2metres_2012).^2));
RMSD_per_5meters_2012_finalsim = sqrt(mean((simulated_finalrho_AVGDP5M-AVGD_per_5metres_2012).^2));

RMSD_2012_finalsim = table(RMSD_per_halfmeter_2012_finalsim,...
    RMSD_per_meter_2012_finalsim,RMSD_per_2meters_2012_finalsim,RMSD_per_5meters_2012_finalsim);
filename_append = append('RMSD_',Year,'_Finsim all 20m sqrt');
writetable(RMSD_2012_finalsim,fullfile(fpath,filename_append));

%% RMSD between final simulated profile and 2012 profile Top 10m at KANU by region

RMSD_TOP5m_2012_finalsim_twentycm_res = sqrt(mean((simulated_finalrho_AVGDP20cm(1:25)-AVGD_2012_Top10m_20cm_RES(1:25)).^2));
RMSD_TOP5m_2012_finalsim_hfm_res = sqrt(mean((simulated_finalrho_AVGDPHM(1:10)-AVGD_2012_Top10m_hfm_RES(1:10)).^2));
RMSD_TOP5m_2012_finalsim_1metre_res = sqrt(mean((simulated_finalrho_AVGDP1M(1:5)-AVGD_2012_Top10m_1m_RES(1:5)).^2));

RMSD_TOP5to10m_2012_finalsim_twentycm_res = sqrt(mean((simulated_finalrho_AVGDP20cm(26:50)-AVGD_2012_Top10m_20cm_RES(26:50)).^2));
RMSD_TOP5to10m_2012_finalsim_hfm_res = sqrt(mean((simulated_finalrho_AVGDPHM(11:20)-AVGD_2012_Top10m_hfm_RES(11:20)).^2));
RMSD_TOP5to10m_2012_finalsim_1metre_res = sqrt(mean((simulated_finalrho_AVGDP1M(6:10)-AVGD_2012_Top10m_1m_RES(6:10)).^2));

RMSD_Top10m_2012_finalsim_twentycm_res = sqrt(mean((simulated_finalrho_AVGDP20cm(1:50)-AVGD_2012_Top10m_20cm_RES(1:50)).^2));
RMSD_Top10m_2012_finalsim_hfm_res = sqrt(mean((simulated_finalrho_AVGDPHM(1:20)-AVGD_2012_Top10m_hfm_RES(1:20)).^2));
RMSD_Top10m_2012_finalsim_1metre_res = sqrt(mean((simulated_finalrho_AVGDP1M(1:10)-AVGD_2012_Top10m_1m_RES(1:10)).^2));

RMSD_2012_regions_20cm_res = table(RMSD_TOP5m_2012_finalsim_twentycm_res, RMSD_TOP5to10m_2012_finalsim_twentycm_res,...
    RMSD_Top10m_2012_finalsim_twentycm_res);
writetable(RMSD_2012_regions_20cm_res,fullfile(fpath,'RMSD_2012_regions_20cm_res sqrt'));

RMSD_2012_regions_hfm_res = table(RMSD_TOP5m_2012_finalsim_hfm_res, RMSD_TOP5to10m_2012_finalsim_hfm_res,...
    RMSD_Top10m_2012_finalsim_hfm_res);
writetable(RMSD_2012_regions_hfm_res,fullfile(fpath,'RMSD_2012_regions_hfm_res sqrt'));

RMSD_2012_regions_1metre_res = table(RMSD_TOP5m_2012_finalsim_1metre_res, RMSD_TOP5to10m_2012_finalsim_1metre_res,...
    RMSD_Top10m_2012_finalsim_1metre_res);
writetable(RMSD_2012_regions_1metre_res,fullfile(fpath,'RMSD_2012_regions_1metre_res sqrt'));

%% RMSD between final simulated profile and 2013 profile at KANU

RMSD_per_halfmeter_2013_finalsim = sqrt(mean((simulated_finalrho_AVGDPHM-AVGD_per_halfmetre_2013).^2));
RMSD_per_meter_2013_finalsim = sqrt(mean((simulated_finalrho_AVGDP1M-AVGD_per_metre_2013).^2));
RMSD_per_2meters_2013_finalsim = sqrt(mean((simulated_finalrho_AVGDP2M-AVGD_per_2metres_2013).^2));
RMSD_per_5meters_2013_finalsim = sqrt(mean((simulated_finalrho_AVGDP5M-AVGD_per_5metres_2013).^2));

RMSD_2013_finalsim = table(RMSD_per_halfmeter_2013_finalsim,...
    RMSD_per_meter_2013_finalsim,RMSD_per_2meters_2013_finalsim,RMSD_per_5meters_2013_finalsim);
writetable(RMSD_2013_finalsim,fullfile(fpath,'RMSD_2013_Finsim sqrt'));

%% RMSD for specific region of firn 2013 profile

fpathRMSD = './toms/PerbsWD89_19/RMSD';

RMSD_TOP5m_2013_finalsim_twentycm_res = sqrt(mean((simulated_finalrho_AVGDP20cm(1:25)-AVGD_per_20cm_Top5m_2013).^2));
RMSD_TOP5m_2013_finalsim_hfm_res = sqrt(mean((simulated_finalrho_AVGDPHM(1:10)-AVGD_per_hfmetre_Top5m_2013).^2));
RMSD_TOP5m_2013_finalsim_1metre_res = sqrt(mean((simulated_finalrho_AVGDP1M(1:5)-AVGD_per_metre_Top5m_2013).^2));

RMSD_TOP5to10m_2013_finalsim_twentycm_res = sqrt(mean((simulated_finalrho_AVGDP20cm(26:50)-AVGD_per_20cm_10to15m_2013).^2));
RMSD_TOP5to10m_2013_finalsim_hfm_res = sqrt(mean((simulated_finalrho_AVGDPHM(11:20)-AVGD_per_hfmetre_10to15m_2013).^2));
RMSD_TOP5to10m_2013_finalsim_1metre_res = sqrt(mean((simulated_finalrho_AVGDP1M(6:10)-AVGD_per_metre_10to15m_2013).^2));

RMSD_TOP10to15m_2013_finalsim_twentycm_res = sqrt(mean((simulated_finalrho_AVGDP20cm(51:75)-AVGD_per_20cm_10to15m_2013).^2));
RMSD_TOP10to15m_2013_finalsim_hfm_res = sqrt(mean((simulated_finalrho_AVGDPHM(21:30)-AVGD_per_hfmetre_10to15m_2013).^2));
RMSD_TOP10to15m_2013_finalsim_1metre_res = sqrt(mean((simulated_finalrho_AVGDP1M(11:15)-AVGD_per_metre_10to15m_2013).^2));

RMSD_TOP15to19m_2013_finalsim_twentycm_res = sqrt(mean((simulated_finalrho_AVGDP20cm(76:95)-AVGD_per_20cm_15to19m_2013)).^2);
RMSD_TOP15to19m_2013_finalsim_hfm_res = sqrt(mean((simulated_finalrho_AVGDPHM(31:38)-AVGD_per_hfmetre_15to19m_2013).^2));
RMSD_15to19m_2013_finalsim_1metre_res = sqrt(mean((simulated_finalrho_AVGDP1M(16:19)-AVGD_per_metre_15to19m_2013).^2));

RMSD_2013_regions_20cm_res = table(RMSD_TOP5m_2013_finalsim_twentycm_res, RMSD_TOP5to10m_2013_finalsim_twentycm_res,...
    RMSD_TOP10to15m_2013_finalsim_twentycm_res,RMSD_TOP15to19m_2013_finalsim_twentycm_res);
writetable(RMSD_2013_regions_20cm_res,fullfile(fpath,'RMSD_2013_regions_20cm_res sqrt'));

RMSD_2013_regions_hfm_res = table(RMSD_TOP5m_2013_finalsim_hfm_res, RMSD_TOP5to10m_2013_finalsim_hfm_res,...
    RMSD_TOP10to15m_2013_finalsim_hfm_res,RMSD_TOP15to19m_2013_finalsim_hfm_res);
writetable(RMSD_2013_regions_hfm_res,fullfile(fpath,'RMSD_2013_regions_hfm_res sqrt'));

RMSD_2013_regions_1metre_res = table(RMSD_TOP5m_2013_finalsim_1metre_res, RMSD_TOP5to10m_2013_finalsim_1metre_res,...
    RMSD_TOP10to15m_2013_finalsim_1metre_res,RMSD_15to19m_2013_finalsim_1metre_res);
writetable(RMSD_2013_regions_1metre_res,fullfile(fpath,'RMSD_2013_regions_1metre_res sqrt'));


%% Plot and data

% initial and final density profiles

figure('DefaultLegendFontSize',Txtsz);
hold on
%s1 = stairs(rho_sim_final(1:100),depth_layers_final(1:100),'Color',[d]);
%s2 = stairs(rho_sim_init(1:100),depth_layers_init(1:100),'Color',[c]);
%s3 = plot(x1_13,y1_13);
%s3 = plot(x1_89,y1_89);
s3 = plot(x1_12,y1_12);
%s3 = plot(x1_12,y1_12);
LWidth =2;
s1.LineWidth = LWidth;
s2.LineWidth = LWidth;
s3.LineWidth = LWidth;
s1.LineStyle = '--';
s2.LineStyle = '--';
%s3.LineStyle = '--';

grid on
ax = gca;
ax.FontSize = Txtsz;
box(ax,'on')
ax.LineWidth = 1.5;

set(ax,'Ydir','reverse')
hold on
%titled = append('Density profiles for ',heads,' perturbation ',Year);
title('Composite Density profile KanU (2012 top 10m) and site J 1989')
%leg_append_sim = append(Year,' Simulated');
%leg_append_meas = append(Year,' Measured');
%legend(leg_append_sim,'1989',leg_append_meas,'Location','west')
xlabel('Density (kg/m^{-3})','Interpreter','tex')
ylabel('Depth(m)')
ylim([0 60])

set(0,'defaultfigurepapersize',[29.7 16]);
set(0,'defaultfigurepaperposition',[.25 .25 [29.7 16]-0.5]);
set(0,'DefaultTextInterpreter','none');
set(0, 'DefaultFigureUnits', 'centimeters');
set(0, 'DefaultFigurePosition', [.25 .25 [29.7 16]-0.5]);
%print(gcf,fullfile(fpath,titled),'-djpeg');
fpath = './toms/PerbsWD89_13';
print(gcf,fullfile(fpath,'Composite Density profile KanU (2012 top 10m and site J 1989'),'-djpeg')
hold off

%% Plot distribution of layers per depth

figure('DefaultLegendFontSize',Txtsz);
hold on
histogram(depth_layers_final,101,'Facecolor',[c]);
histogram(depth_layers_init,101,'Facecolor',[d]);
T20layers_fin = depth_layers_final<=20;
T20layers_init = depth_layers_init <=20;
T20Layerscount_fin = sum(T20layers_fin);
T20Layerscount_init = sum(T20layers_init);

grid on
ax = gca;

box(ax,'on')
ax.LineWidth = 1.5;
ax.FontSize = Txtsz;
titled_layer = append('Distribution of layers, top 20m of firn for ',heads,' perturbation');
title(titled_layer)
xlabel('Depth (m)')
ylabel('Count')

legend(Year, '1989')
str_append = append('Number of layers in ',Year,' =');
str = {str_append,T20Layerscount_fin, 'Number of layers in 1989 =',T20Layerscount_init};
text(20,5,str,'FontSize',Txtsz)

set(0,'defaultfigurepapersize',[29.7 16]);
set(0,'defaultfigurepaperposition',[.25 .25 [29.7 16]-0.5]);
set(0,'DefaultTextInterpreter','none');
set(0, 'DefaultFigureUnits', 'centimeters');
set(0, 'DefaultFigurePosition', [.25 .25 [29.7 16]-0.5]);
print(gcf,fullfile(fpath,titled_layer),'-djpeg');

hold off

%% Interpolation plot checker final time step sim

figure('DefaultLegendFontSize',Txtsz);
fig = gcf;

subplot(1,2,1)
I1 = stairs(Top20_Densities_final,Top20_Depth_final,'Color',[c]);
 
grid on

I1.LineWidth = LW;
ax = gca;
box(ax,'on')
ax.LineWidth = 1.5;
ax.FontSize = Txtsz;

set(ax,'Ydir','reverse')
title_appended = append ('Simulated profile at final time-step ',Year);
title(title_appended,'FontSize',Txtsz)
xlabel('Density (kg m^{⁻3})','Interpreter','tex')
ylabel('Depth(m)')
ylim([0 20])
subplot(1,2,2)
plot(y_Densities_Final,x_Depths_Final,'o',yi_interpolated_densities_final,xi_interpolated_depths_final,...
    'MarkerFaceColor',[d],'MarkerEdgeColor','k','Color',[c]);
set(gca,'Ydir','reverse')

grid on

ax = gca;
box(ax,'on')
ax.LineWidth = 1.5;
ax.FontSize = Txtsz;
title('Interpolated profile')
xlabel('Density (kg m^{⁻3})','Interpreter','tex')
ylabel('Depth (m)')
ylim([0 20])
Leg = legend('Interpolated points','Location','west');
set(Leg,'Color','none')
set(0,'defaultfigurepapersize',[29.7 16]);
set(0,'defaultfigurepaperposition',[.25 .25 [29.7 16]-0.5]);
set(0,'DefaultTextInterpreter','none');
set(0, 'DefaultFigureUnits', 'centimeters');
set(0, 'DefaultFigurePosition', [.25 .25 [29.7 16]-0.5]);
filename = append('Simulated profile & interpolated data ',Year);
print(gcf,fullfile(fpath,filename),'-djpeg');


set(findall(gcf,'-property','FontSize'),'FontSize',Txtsz)
saveas(gcf,fullfile(fpath,'Simulated density profile, Top 50m, final timestep & interpolated data'),'pdf');

% Interpolation plot checker initial time step sim


figure('DefaultLegendFontSize',Txtsz);
subplot(1,2,1)
II = stairs(Top20_Densities_init,Top20_Depth_init,'Color',[c]);

grid on
II.LineWidth = LW;
ax = gca;
box(ax,'on')
ax.LineWidth = 1.5;
ax.FontSize = Txtsz;
set(ax,'Ydir','reverse')
title('Simulated profile at initial time-step 1989')
xlabel('Density (kg m^{⁻3})','Interpreter','tex')
ylabel('Depth(m)')
ylim([0 20])
subplot(1,2,2)
plot(y_Densities_init,x_Depths_init,'o',yi_interpolated_densities_init,xi_interpolated_depths_init,...
    'MarkerFaceColor',[d],'MarkerEdgeColor','k','Color',[c]);

grid on
ax = gca;
box(ax,'on')
ax.LineWidth = 1.5;
ax.FontSize = Txtsz;
set(ax,'Ydir','reverse')
title('Interpolated profile')
xlabel('Density (kg m^{⁻3})','Interpreter','tex')
ylabel('Depth (m)')
ylim([0 20])
Leg = legend('Interpolated data points','Location','west');
set(Leg,'Color','none')
set(0,'defaultfigurepapersize',[29.7 16]);
set(0,'defaultfigurepaperposition',[.25 .25 [29.7 16]-0.5]);
set(0,'DefaultTextInterpreter','none');
set(0, 'DefaultFigureUnits', 'centimeters');
set(0, 'DefaultFigurePosition', [.25 .25 [29.7 16]-0.5]);
filename_appended = append('Simulated profile at initial timestep 1989 & interpolated data ',(Year));
print(gcf,fullfile(fpath,filename_appended),'-djpeg');

%% Plot site J against perturbation

DP_1989 = readtable('DensityProfile_KAN-U_1989_metres.csv');
DP89_depths = table2array(DP_1989(:,1));
DP89_densities = table2array(DP_1989(:,2));

compare with perturbed 
Pert = append(pert,number,'.csv');
DP_1989_Pert = readtable(Pert);
DP89_depths_Pert = table2array(DP_1989_Pert(:,1));
DP89_densities_Pert = table2array(DP_1989_Pert(:,2));

figure('DefaultLegendFontSize',Txtsz);
D1 = stairs(DP89_densities,DP89_depths);
D1.LineWidth = LW;

grid on

ax = gca;
box(ax,'on')
ax.LineWidth = 1.5;
ax.FontSize = Txtsz;
set(ax,'Ydir','reverse')
titled = append('Top 20m at site J with ',Type,number,' perturbation ',Year);
title(titled)
xlabel('Density (kg m^{⁻3})','Interpreter','tex')
ylabel('Depth (m)')
cm ice slab added to each metre of firn
ylim([0 10])

hold on 
D2 = stairs(DP89_densities_Pert,DP89_depths_Pert,'--','Color',[e]);
D2.LineWidth = LW;
legend_append = append(Type,number,' perturbation');
legend('Site J', legend_append)
set(0,'defaultfigurepapersize',[29.7 16]);
set(0,'defaultfigurepaperposition',[.25 .25 [29.7 16]-0.5]);
set(0,'DefaultTextInterpreter','none');
set(0, 'DefaultFigureUnits', 'centimeters');
set(0, 'DefaultFigurePosition', [.25 .25 [29.7 16]-0.5]);
filename_append = append('Site J and perturbed profile ', heads, Year);
print(gcf,fullfile(fpath,filename_append),'-djpeg');

hold off
close all