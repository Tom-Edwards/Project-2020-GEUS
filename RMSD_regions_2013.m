%% RMSD region plotter
% Homo goes from 6:2:20 + 25, 30,35
% Thick goes through 10:2:30

%RMSD_2013_main = readtable('RMSD_.csv');

number = '30';      % perturbation number
pert = 'Thick';     % Thick or Homo, (slab addition or homogeneous perturbation)
heads = append(pert,number);
outputFP = append('/toms/PerbsWD89_13/P',heads);

fpath = './toms/PerbsWD89_13/RMSD/Regions_2013';

twentycm_res = append(outputFP,'/RMSD_2013_regions_20cm_res sqrt.txt');
halfmetre_res = append(outputFP,'/RMSD_2013_regions_hfm_res sqrt.txt');
onemetre_res = append(outputFP,'/RMSD_2013_regions_1metre_res sqrt.txt');

%% Ice slab variables
% 
% Thick10_20cm = table2array(readtable(twentycm_res));
% Thick10_hfm_res = table2array(readtable(halfmetre_res));
% Thick10_1m = table2array(readtable(onemetre_res));
 
% Thick12_20cm = table2array(readtable(twentycm_res));
% Thick12_hfm_res = table2array(readtable(halfmetre_res));
% Thick12_1m = table2array(readtable(onemetre_res));

% Thick14_20cm = table2array(readtable(twentycm_res));
% Thick14_hfm_res = table2array(readtable(halfmetre_res));
% Thick14_1m = table2array(readtable(onemetre_res));
% 
% Thick16_20cm = table2array(readtable(twentycm_res));
% Thick16_hfm_res = table2array(readtable(halfmetre_res));
% Thick16_1m = table2array(readtable(onemetre_res));

% Thick18_20cm = table2array(readtable(twentycm_res));
% Thick18_hfm_res = table2array(readtable(halfmetre_res));
% Thick18_1m = table2array(readtable(onemetre_res));
% 
% Thick20_20cm = table2array(readtable(twentycm_res));
% Thick20_hfm_res = table2array(readtable(halfmetre_res));
% Thick20_1m = table2array(readtable(onemetre_res));

% Thick22_20cm = table2array(readtable(twentycm_res));
% Thick22_hfm_res = table2array(readtable(halfmetre_res));
% Thick22_1m = table2array(readtable(onemetre_res));
% 
% Thick24_20cm = table2array(readtable(twentycm_res));
% Thick24_hfm_res = table2array(readtable(halfmetre_res));
% Thick24_1m = table2array(readtable(onemetre_res));
% 
% Thick26_20cm = table2array(readtable(twentycm_res));
% Thick26_hfm_res = table2array(readtable(halfmetre_res));
% Thick26_1m = table2array(readtable(onemetre_res));

% Thick28_20cm = table2array(readtable(twentycm_res));
% Thick28_hfm_res = table2array(readtable(halfmetre_res));
% Thick28_1m = table2array(readtable(onemetre_res));

% Thick30_20cm = table2array(readtable(twentycm_res));
% Thick30_hfm_res = table2array(readtable(halfmetre_res));
% Thick30_1m = table2array(readtable(onemetre_res));

%% Homogeneous Variables
% 
% Homo6_20cm = table2array(readtable(twentycm_res));
% Homo6_hfmetre_res = table2array(readtable(halfmetre_res));
% Homo6_1m = table2array(readtable(onemetre_res));

% Homo8_20cm = table2array(readtable(twentycm_res));
% Homo8_hfmetre_res = table2array(readtable(halfmetre_res));
% Homo8_1m = table2array(readtable(onemetre_res));
% 
% Homo10_20cm = table2array(readtable(twentycm_res));
% Homo10_hfmetre_res = table2array(readtable(halfmetre_res));
% Homo10_1m = table2array(readtable(onemetre_res));
% 
% Homo12_20cm = table2array(readtable(twentycm_res));
% Homo12_hfmetre_res = table2array(readtable(halfmetre_res));
% Homo12_1m = table2array(readtable(onemetre_res));
% 
% Homo14_20cm = table2array(readtable(twentycm_res));
% Homo14_hfmetre_res = table2array(readtable(halfmetre_res));
% Homo14_1m = table2array(readtable(onemetre_res));
% 
% Homo16_20cm = table2array(readtable(twentycm_res));
% Homo16_hfmetre_res = table2array(readtable(halfmetre_res));
% Homo16_1m = table2array(readtable(onemetre_res));
% 
% Homo18_20cm = table2array(readtable(twentycm_res));
% Homo18_hfmetre_res = table2array(readtable(halfmetre_res));
% Homo18_1m = table2array(readtable(onemetre_res));
% 
% Homo20_20cm = table2array(readtable(twentycm_res));
% Homo20_hfmetre_res = table2array(readtable(halfmetre_res));
% Homo20_1m = table2array(readtable(onemetre_res));

% Homo25_20cm = table2array(readtable(twentycm_res));
% Homo25_hfmetre_res = table2array(readtable(halfmetre_res));
% Homo25_1m = table2array(readtable(onemetre_res));
% 
% Homo30_20cm = table2array(readtable(twentycm_res));
% Homo30_hfmetre_res = table2array(readtable(halfmetre_res));
% Homo30_1m = table2array(readtable(onemetre_res));

% Homo35_20cm = table2array(readtable(twentycm_res));
% Homo35_hfmetre_res = table2array(readtable(halfmetre_res));
% Homo35_1m = table2array(readtable(onemetre_res));

%% Homogeneous perturbations ---------------------------------------------

%% 20cm resolution by region depth

Homo_Fivemetre_region_20cm_res = table(Homo6_20cm(1),Homo8_20cm(1),Homo10_20cm(1),Homo12_20cm(1),Homo14_20cm(1),Homo16_20cm(1),...
    Homo18_20cm(1),Homo20_20cm(1),Homo25_20cm(1),Homo30_20cm(1),Homo35_20cm(1));

Homo_Fivetoten_metre_region_20cm_res = table(Homo6_20cm(2),Homo8_20cm(2),Homo10_20cm(2),Homo12_20cm(2),Homo14_20cm(2),Homo16_20cm(2),...
    Homo18_20cm(2),Homo20_20cm(2),Homo25_20cm(2),Homo30_20cm(2),Homo35_20cm(2));

Homo_Tentofifteen_metre_region_20cm_res = table(Homo6_20cm(3),Homo8_20cm(3),Homo10_20cm(3),Homo12_20cm(3),Homo14_20cm(3),Homo16_20cm(3),...
    Homo18_20cm(3),Homo20_20cm(3),Homo25_20cm(3),Homo30_20cm(3),Homo35_20cm(3));

Homo_Fifteentonineteen_metre_region_20cm_res = table(Homo6_20cm(4),Homo8_20cm(1),Homo10_20cm(4),Homo12_20cm(4),Homo14_20cm(4),Homo16_20cm(4),...
    Homo18_20cm(4),Homo20_20cm(4),Homo25_20cm(4),Homo30_20cm(4),Homo35_20cm(4));

Homo_Fivemetre_region_20cm_res = table2array(Homo_Fivemetre_region_20cm_res)';
Homo_Fivetoten_metre_region_20cm_res = table2array(Homo_Fivetoten_metre_region_20cm_res)';
Homo_Tentofifteen_metre_region_20cm_res = table2array(Homo_Tentofifteen_metre_region_20cm_res)';
Homo_Fifteentonineteen_metre_region_20cm_res = table2array(Homo_Fifteentonineteen_metre_region_20cm_res)';


Twenty_cm_res_regions_Homo = table(Homo_Fivemetre_region_20cm_res,Homo_Fivetoten_metre_region_20cm_res,Homo_Tentofifteen_metre_region_20cm_res,...
    Homo_Fifteentonineteen_metre_region_20cm_res);

writetable(Twenty_cm_res_regions_Homo,fullfile(fpath,'20cm resolution, homo regionals 2013 sqrt'));

%% 50cm resoution by region depth

Homo_Fivemetre_region_halfmetre_res = table(Homo6_hfmetre_res(1),Homo8_hfmetre_res(1),Homo10_hfmetre_res(1),Homo12_hfmetre_res(1),Homo14_hfmetre_res(1),Homo16_hfmetre_res(1),...
    Homo18_hfmetre_res(1),Homo20_hfmetre_res(1),Homo25_hfmetre_res(1),Homo30_hfmetre_res(1),Homo35_hfmetre_res(1));

Homo_Fivetoten_metre_region_halfmetre_res = table(Homo6_hfmetre_res(2),Homo8_hfmetre_res(2),Homo10_hfmetre_res(2),Homo12_hfmetre_res(2),Homo14_hfmetre_res(2),Homo16_hfmetre_res(2),...
    Homo18_hfmetre_res(2),Homo20_hfmetre_res(2),Homo25_hfmetre_res(2),Homo30_hfmetre_res(2),Homo35_hfmetre_res(2));

Homo_Tentofifteen_metre_region_halfmetre_res = table(Homo6_hfmetre_res(3),Homo8_hfmetre_res(3),Homo10_hfmetre_res(3),Homo12_hfmetre_res(3),Homo14_hfmetre_res(3),Homo16_hfmetre_res(3),...
    Homo18_hfmetre_res(3),Homo20_hfmetre_res(3),Homo25_hfmetre_res(3),Homo30_hfmetre_res(3),Homo35_hfmetre_res(3));

Homo_Fifteentonineteen_metre_region_halfmetre_res = table(Homo6_hfmetre_res(4),Homo8_hfmetre_res(1),Homo10_hfmetre_res(4),Homo12_hfmetre_res(4),Homo14_hfmetre_res(4),Homo16_hfmetre_res(4),...
    Homo18_hfmetre_res(4),Homo20_hfmetre_res(4),Homo25_hfmetre_res(4),Homo30_hfmetre_res(4),Homo35_hfmetre_res(4));

Homo_Fivemetre_region_halfmetre_res = table2array(Homo_Fivemetre_region_halfmetre_res)';
Homo_Fivetoten_metre_region_halfmetre_res = table2array(Homo_Fivetoten_metre_region_halfmetre_res)';
Homo_Tentofifteen_metre_region_halfmetre_res = table2array(Homo_Tentofifteen_metre_region_halfmetre_res)';
Homo_Fifteentonineteen_metre_region_halfmetre_res = table2array(Homo_Fifteentonineteen_metre_region_halfmetre_res)';

Fifty_cm_res_regions_Homo = table(Homo_Fivemetre_region_halfmetre_res,Homo_Fivetoten_metre_region_halfmetre_res,Homo_Tentofifteen_metre_region_halfmetre_res,...
    Homo_Fifteentonineteen_metre_region_halfmetre_res);

writetable(Fifty_cm_res_regions_Homo,fullfile(fpath,'50cm resolution, homo regionals 2013 sqrt'));


%% 1m resoution by region depth

Homo_Fivemetre_region_metre_res = table(Homo6_1m(1),Homo8_1m(1),Homo10_1m(1),Homo12_1m(1),Homo14_1m(1),Homo16_1m(1),...
    Homo18_1m(1),Homo20_1m(1),Homo25_1m(1),Homo30_1m(1),Homo35_1m(1));

Homo_Fivetoten_metre_region_metre_res = table(Homo6_1m(2),Homo8_1m(2),Homo10_1m(2),Homo12_1m(2),Homo14_1m(2),Homo16_1m(2),...
    Homo18_1m(2),Homo20_1m(2),Homo25_1m(2),Homo30_1m(2),Homo35_1m(2));

Homo_Tentofifteen_metre_region_metre_res = table(Homo6_1m(3),Homo8_1m(3),Homo10_1m(3),Homo12_1m(3),Homo14_1m(3),Homo16_1m(3),...
    Homo18_1m(3),Homo20_1m(3),Homo25_1m(3),Homo30_1m(3),Homo35_1m(3));

Homo_Fifteentonineteen_metre_region_metre_res = table(Homo6_1m(4),Homo8_1m(1),Homo10_1m(4),Homo12_1m(4),Homo14_1m(4),Homo16_1m(4),...
    Homo18_1m(4),Homo20_1m(4),Homo25_1m(4),Homo30_1m(4),Homo35_1m(4));

Homo_Fivemetre_region_metre_res = table2array(Homo_Fivemetre_region_metre_res)';
Homo_Fivetoten_metre_region_metre_res = table2array(Homo_Fivetoten_metre_region_metre_res)';
Homo_Tentofifteen_metre_region_metre_res = table2array(Homo_Tentofifteen_metre_region_metre_res)';
Homo_Fifteentonineteen_metre_region_metre_res = table2array(Homo_Fifteentonineteen_metre_region_metre_res)';


One_metre_res_regions_Homo = table(Homo_Fivemetre_region_metre_res,Homo_Fivetoten_metre_region_metre_res,Homo_Tentofifteen_metre_region_metre_res,...
    Homo_Fifteentonineteen_metre_region_metre_res);

writetable(One_metre_res_regions_Homo,fullfile(fpath,'1m resolution, homo regionals 2013 sqrt'));


%% Slab perturbations ---------------------------------------------

%% 20cm resolution by region depth

Thick_Fivemetre_region_20cm_res = table(Thick10_20cm(1),Thick12_20cm(1),Thick14_20cm(1),Thick16_20cm(1),...
    Thick18_20cm(1),Thick20_20cm(1),Thick22_20cm(1),Thick24_20cm(1),Thick26_20cm(1),Thick28_20cm(1),Thick30_20cm(1));

Thick_Fivetoten_metre_region_20cm_res = table(Thick10_20cm(2),Thick12_20cm(2),Thick14_20cm(2),Thick16_20cm(2),...
    Thick18_20cm(2),Thick20_20cm(2),Thick22_20cm(2),Thick24_20cm(2),Thick26_20cm(2),Thick28_20cm(2),Thick30_20cm(2));

Thick_Tentofifteen_metre_region_20cm_res = table(Thick10_20cm(3),Thick12_20cm(3),Thick14_20cm(3),Thick16_20cm(3),...
    Thick18_20cm(3),Thick20_20cm(3),Thick22_20cm(3),Thick24_20cm(3),Thick26_20cm(3),Thick28_20cm(3),Thick30_20cm(3));

Thick_Fifteentonineteen_metre_region_20cm_res = table(Thick10_20cm(4),Thick12_20cm(4),Thick14_20cm(4),Thick16_20cm(4),...
    Thick18_20cm(4),Thick20_20cm(4),Thick22_20cm(4),Thick24_20cm(4),Thick26_20cm(4),Thick28_20cm(4),Thick30_20cm(4));

Thick_Fivemetre_region_20cm_res = table2array(Thick_Fivemetre_region_20cm_res)';
Thick_Fivetoten_metre_region_20cm_res = table2array(Thick_Fivetoten_metre_region_20cm_res)';
Thick_Tentofifteen_metre_region_20cm_res = table2array(Thick_Tentofifteen_metre_region_20cm_res)';
Thick_Fifteentonineteen_metre_region_20cm_res = table2array(Thick_Fifteentonineteen_metre_region_20cm_res)';


Twenty_cm_res_regions_Thick = table(Thick_Fivemetre_region_20cm_res,Thick_Fivetoten_metre_region_20cm_res,Thick_Tentofifteen_metre_region_20cm_res,...
    Thick_Fifteentonineteen_metre_region_20cm_res);

writetable(Twenty_cm_res_regions_Thick,fullfile(fpath,'20cm resolution, thick regionals 2013 sqrt'));



%% 50cm resoution by region depth

Thick_Fivemetre_region_halfmetre_res = table(Thick10_hfm_res(1),Thick12_hfm_res(1),Thick14_hfm_res(1),Thick16_hfm_res(1),...
    Thick18_hfm_res(1),Thick20_hfm_res(1),Thick22_hfm_res(1),Thick24_hfm_res(1),Thick26_hfm_res(1),Thick28_hfm_res(1),Thick30_hfm_res(1));

Thick_Fivetoten_metre_region_halfmetre_res = table(Thick10_hfm_res(2),Thick12_hfm_res(2),Thick14_hfm_res(2),Thick16_hfm_res(2),...
    Thick18_hfm_res(2),Thick20_hfm_res(2),Thick22_hfm_res(2),Thick24_hfm_res(2),Thick26_hfm_res(2),Thick28_hfm_res(2),Thick30_hfm_res(2));

Thick_Tentofifteen_metre_region_halfmetre_res = table(Thick10_hfm_res(3),Thick12_hfm_res(3),Thick14_hfm_res(3),Thick16_hfm_res(3),...
    Thick18_hfm_res(3),Thick20_hfm_res(3),Thick22_hfm_res(3),Thick24_hfm_res(3),Thick26_hfm_res(3),Thick28_hfm_res(3),Thick30_hfm_res(3));

Thick_Fifteentonineteen_metre_region_halfmetre_res = table(Thick10_hfm_res(4),Thick12_hfm_res(1),Thick14_hfm_res(4),Thick16_hfm_res(4),...
    Thick18_hfm_res(4),Thick20_hfm_res(4),Thick22_hfm_res(4),Thick24_hfm_res(4),Thick26_hfm_res(4),Thick28_hfm_res(4),Thick30_hfm_res(4));

Thick_Fivemetre_region_halfmetre_res = table2array(Thick_Fivemetre_region_halfmetre_res)';
Thick_Fivetoten_metre_region_halfmetre_res = table2array(Thick_Fivetoten_metre_region_halfmetre_res)';
Thick_Tentofifteen_metre_region_halfmetre_res = table2array(Thick_Tentofifteen_metre_region_halfmetre_res)';
Thick_Fifteentonineteen_metre_region_halfmetre_res = table2array(Thick_Fifteentonineteen_metre_region_halfmetre_res)';

Halfmetre_cm_res_regions_Thick = table(Thick_Fivemetre_region_halfmetre_res,Thick_Fivetoten_metre_region_halfmetre_res,Thick_Tentofifteen_metre_region_halfmetre_res,...
   Thick_Fifteentonineteen_metre_region_halfmetre_res);

writetable(Halfmetre_cm_res_regions_Thick,fullfile(fpath,'50cm resolution, thick regionals 2013 sqrt'));



%% 1m resoution by region depth

Thick_Fivemetre_region_metre_res = table(Thick10_1m(1),Thick12_1m(1),Thick14_1m(1),Thick16_1m(1),...
    Thick18_1m(1),Thick20_1m(1),Thick22_1m(1),Thick24_1m(1),Thick26_1m(1),Thick28_1m(1),Thick30_1m(1));

Thick_Fivetoten_metre_region_metre_res = table(Thick10_1m(2),Thick12_1m(2),Thick14_1m(2),Thick16_1m(2),...
    Thick18_1m(2),Thick20_1m(2),Thick22_1m(2),Thick24_1m(2),Thick26_1m(2),Thick28_1m(2),Thick30_1m(2));

Thick_Tentofifteen_metre_region_metre_res = table(Thick10_1m(3),Thick12_1m(3),Thick14_1m(3),Thick16_1m(3),...
    Thick18_1m(3),Thick20_1m(3),Thick22_1m(3),Thick24_1m(3),Thick26_1m(3),Thick28_1m(3),Thick30_1m(3));

Thick_Fifteentonineteen_metre_region_metre_res = table(Thick10_1m(4),Thick12_1m(4),Thick14_1m(4),Thick16_1m(4),...
    Thick18_1m(4),Thick20_1m(4),Thick22_1m(4),Thick24_1m(4),Thick26_1m(4),Thick28_1m(4),Thick30_1m(4));

Thick_Fivemetre_region_metre_res = table2array(Thick_Fivemetre_region_metre_res)';
Thick_Fivetoten_metre_region_metre_res = table2array(Thick_Fivetoten_metre_region_metre_res)';
Thick_Tentofifteen_metre_region_metre_res = table2array(Thick_Tentofifteen_metre_region_metre_res)';
Thick_Fifteentonineteen_metre_region_metre_res = table2array(Thick_Fifteentonineteen_metre_region_metre_res)';

Metre_cm_res_regions_Thick = table(Thick_Fivemetre_region_metre_res,Thick_Fivetoten_metre_region_metre_res,Thick_Tentofifteen_metre_region_metre_res,...
   Thick_Fifteentonineteen_metre_region_metre_res);

writetable(Metre_cm_res_regions_Thick,fullfile(fpath,'1m resolution, thick regionals 2013 sqrt'));


%% Plotting 
%import files

Twenties = table2array(Twenty_cm_res_regions_Homo);
Fifties = table2array(Fifty_cm_res_regions_Homo);
Metries = table2array(One_metre_res_regions_Homo);

Twenties_Thick = table2array(Twenty_cm_res_regions_Thick);
Fifties_Thick = table2array(Halfmetre_cm_res_regions_Thick);
Metries_Thick = table2array(Metre_cm_res_regions_Thick);
%%

thick_index = 10:2:30;
homo_index = 6:2:20;
homo_index = [homo_index 25 30 35];

ind1 = homo_index;
ind2 = thick_index;

%%

figure;
grid on
hold on
Year = ' 2013';

res = '20cm';
type = 'Thick';

s1 = plot(ind2,Twenties_Thick(:,1));
s2 = plot(ind2,Twenties_Thick(:,2));
s3 = plot(ind2,Twenties_Thick(:,3));
% s1 = plot(ind2,Fifties_Thick(:,1));
% s2 = plot(ind2,Fifties_Thick(:,2));
% s3 = plot(ind2,Fifties_Thick(:,3));
% s1 = plot(ind2,Metries_Thick(:,1));
% s2 = plot(ind2,Metries_Thick(:,2));
% s3 = plot(ind2,Metries_Thick(:,3));
% 
% s1 = plot(ind1,Twenties(:,1));
% s2 = plot(ind1,Twenties(:,2));
% s3 = plot(ind1,Twenties(:,3));
% s1 = plot(ind1,Fifties(:,1));
% s2 = plot(ind1,Fifties(:,2));
% s3 = plot(ind1,Fifties(:,3));
% s1 = plot(ind1,Metries(:,1));
% s2 = plot(ind1,Metries(:,2));
% s3 = plot(ind1,Metries(:,3));

Txtsz = 16;
LW = 2; % Line width for all plots
ax = gca;
box(ax,'on')
ax.LineWidth = LW;
ax.FontSize = Txtsz;
s1.LineWidth = LW;
s2.LineWidth = LW;
s3.LineWidth = LW;
% 
% xticks(homo_index)
% xticklabels({'6','8','10','12','14','16','18','20','25','30','35'})
% xlabel('Homogenous perturbation(%)')
% ylabel('RMSD 2013 - simulation (kg/m^{-3})','Interpreter','tex')
% title_appended = append('RMSD ',type,' perturbations ',res,Year);
% title(title_appended)
% HLegend = legend({'0:5m','5:10m','10:15m'},'Location','northeast');
% set(HLegend,'Color','none')
% %ylim([90 330])

xticks(thick_index)
xticklabels({'10','12','14','16','18','20','22','24','26','28','30'})
xlabel('Thick perturbation(cm)')
ylabel('RMSD 2013 - simulation (kg/m^{-3})','Interpreter','tex')
title_appended = append('RMSD ',type,' perturbations ',res,Year);
title(title_appended)
HLegend = legend({'0:5m','5:10m','10:15m',},'Location','northeast');
set(HLegend,'Color','none')
%ylim([90 340])

set(0,'defaultfigurepapersize',[29.7 16]);
set(0,'defaultfigurepaperposition',[.25 .25 [29.7 16]-0.5]);
set(0,'DefaultTextInterpreter','none');
set(0, 'DefaultFigureUnits', 'centimeters');
set(0, 'DefaultFigurePosition', [.25 .25 [29.7 16]-0.5]);

filename_appended = append('RMSD ',type,' perturbations ',res,Year);

print(gcf,fullfile(fpath,filename_appended),'-djpeg')
hold off