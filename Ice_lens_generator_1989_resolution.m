%% Import profile
clear
clc
DP_1989 = readtable('DensityProfile_KAN-U_1989_metres.csv');

%% Random ice slab generator in the upper 10 metres of firn

% Take top 10m subtract difference between density at each level and max
% density

Base = DP_1989{1:1200,2};
Diff = zeros(1000,1);
for i = 1:length(Base)
    Diff(i) = 917 - Base(i);
end

% Pick a random number that generates a random number of ice lenses
rl = round(100.*rand(1,1));  
% Index positions for each lens (number between 10:1000)
index = round((1000-10).*rand(rl,1) + 10);
% empty array for new density values
new_density = zeros(1200,1);                
% Generate new density values 
for i = 1:rl
    index(i)
    new_density(index(i)) = (917-Base(index(i))).*rand(1,1) + Base(index(i)) ;
    
end
% Store them somewhere
for i = 1:rl
    Base(index(i),2) = new_density(index(i));  
end

%Convert original density profile from table to an array
DPtab = table2array(DP_1989(:,(1:2)));

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
ylim([0 100])
xlabel('Density(kg m^{⁻3})');
ylabel('Depth(m)');
title('Density profile Site J 1989')
set(gca, 'YDir','reverse') % Reverses direction of axis so surface is top of y axis

%% Export new profile as csv

writetable(DP_1989,'DensityProfile_KAN-U_1989_toospiky2.csv')
 