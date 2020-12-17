%% This script was used to calculate the average density of firn in the top 20m of the KANU 2012 Density profile
% Considerations to be made - Take out the first 0.5 - 1m before makeing the average density calculations. 
% Rewrite script so that average density increaser is in a separate .m file
% written by Tom Edwards 11/2020 t.ron.edwards@gmail.com

%% Take top 20m of 89 profile

DP_1989 = readtable('DensityProfile_KAN-U_1989_metres.csv');
DP89_T20_depths = table2array(DP_1989(1:2000,1));
DP89_T20_densities = table2array(DP_1989(1:2000,2));



%% Plot top20m of 2012 profile

Top20_densities = table2array(DP_2012(1:200,2));
Top20_depths = table2array(DP_2012(1:200,1));


plot(Top20_densities,Top20_depths)
xlabel('Density(kg m^{⁻3})');
ylabel('Depth(m)');
title('KANU 2012 profile top 20m')
set(gca, 'YDir','reverse')

%% Calculate and display average density per 5 measures (pseudo 1/2 metre)

A1 = reshape(Top20,5,40);                       % Top 20m divided into 1/2 metre columns (5 data points)
C1 = sum(A1);                                   % Sum of each column
AVGDhfm_2012 = (C1./5)';                              % Average of every 5 readings
% plot(B1avg);
% B1 = reshape(sum(reshape(Top20,5,40)),[],1);  % 

%% Calculate and display average density per metre 

A = reshape(DP89_T20_densities,10,20);
C = sum(A);
AVGD1m_2012 = (C./10)';
% plot(Bavg)
% B = reshape(sum(reshape(Top20,10,20)),[],1);

%% Calculate and display average density per 2 metres

A2 = reshape(Top20,20,10);
C2 = sum(A2);
AVGD2m_2012 = (C2./20)';
% plot(B2avg)
% B2 = reshape(sum(reshape(Top20,20,10)),[],1);


%% Change average density by some percentage linearly

Newaverage = AVGD1*1.8;
figure
plot(B2avg)
hold on
plot(Newaverage)

%% Redistribute randomly 

m = round(sum(Newaverage)); % Choose the desired sum
n = 20; % Choose the number of integers
v = diff([0,sort(randperm(m-1,n-1)),m]);

plot(v)
repasteme = v';

% v is a vector that can be repasted in to replace the corresponding
% section of firn who's density was changed.

%% Write 2012 average density per metre(s) profile to csv. 

writematrix(AVGD1m_2012, '2012 average density per metre');
writematrix(AVGDhfm_2012, '2012 average density per halfmetre');
writematrix(AVGD2m_2012, '2012 average density per 2metres');

