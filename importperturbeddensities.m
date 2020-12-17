
P0 = readtable('Initial_profile.csv');
P1 = readtable('Homo10.csv');
P2 = readtable('Homo12.csv');
P3 = readtable('Homo14.csv');
P4 = readtable('Homo16.csv');
P5 = readtable('Homo18.csv');
P6 = readtable('Homo20.csv');
P7 = readtable('Homo25.csv');
P8 = readtable('Homo30.csv');
% P7 = readtable('P7.csv');
% P8 = readtable('P8.csv');
% P9 = readtable('P9.csv');
% P10 = readtable('P10.csv');

P0dens = table2array(P0(:,2));

depths = table2array(P1(:,1));
P1dens = table2array(P1(:,2));
P2dens = table2array(P2(:,2));
P3dens = table2array(P3(:,2));
P4dens = table2array(P4(:,2));
P5dens = table2array(P5(:,2));
P6dens = table2array(P6(:,2));
P7dens = table2array(P7(:,2));
P8dens = table2array(P8(:,2));
% P9dens = table2array(P9(:,2));
% P10dens = table2array(P10(:,2));

figure
set(gca,'ydir','reverse')
hold on
%plot(P0dens,depths)
plot(P1dens,depths)
plot(P2dens,depths)
plot(P3dens,depths)
plot(P4dens,depths)
plot(P5dens,depths)
plot(P6dens,depths)
plot(P7dens,depths)
plot(P8dens,depths)
% plot(P9dens,depths)
% plot(P10dens,depths)