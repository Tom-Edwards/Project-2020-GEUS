%% Import profile
clear
clc
DP_2012 = readtable('RetMIP_density_KAN-U 2012.csv');

%% Random ice slab generator in the upper 10 metres of firn

% Take top 10m subtract difference between density at each level and max
% density. Column 2 is densities Column 1 Depths in 10cm increments

Base = DP_2012{1:100,2};
Diff = zeros(100,1);
MaxD = max(DP_2012{1:100,2});
for i = 1:length(Base)
    Diff(i) = MaxD - Base(i);
end

% Pick a random number to generate a random number of ice lenses
rl = round(100.*rand(1,1));  
% Index positions for each lens (number between 10:1000)
index = round((100-10).*rand(rl,1) + 10);
% empty array for new density values
new_density = zeros(100,1);                
% Generate new density values 
for i = 1:rl
    index(i)
    new_density(index(i)) = ((MaxD-Base(index(i))).*rand(1,1)) + Base(index(i)) ;
    
end
% Store them somewhere
for i = 1:rl
    Base(index(i),2) = new_density(index(i));  
end

%Convert original density profile from table to an array
DPtab = table2array(DP_2012(:,(1:2)));

%Insert new density values
for i = 1:rl
    DPtab(index(i),2) = Base(index(i),2);
    
end

% Convert back to table, plot
DP_new = array2table(DPtab);
figure;
x2 = DP_new{:,2};
y2 = DP_new{:,1};
plot(x2,y2)
%ylim([0 100])
xlabel('Density(kg m^{⁻3})');
ylabel('Depth(m)');
title('Density profile KANU 2012')
set(gca, 'YDir','reverse') % Reverses direction of axis so surface is top of y axis

%% Export new profile as csv

writetable(DP_2012,'DensityProfile_KAN-U_1989_toospiky2.csv')
 