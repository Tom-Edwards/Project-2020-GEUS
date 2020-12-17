%% Import profile

DP_2012 = readtable('RetMIP_density_KAN-U 2012.csv');
figure;
x2 = table2array(DP_2012(:,2));
y2 = table2array(DP_2012(:,1));
plot(x2,y2)
%ylim([0 100])
xlabel('Density(kg m^{⁻3})');
ylabel('Depth(m)');
title('Random density profile')
set(gca, 'YDir','reverse')

%% Pick 150 random densities from the 2012 profile

DP2012 = table2array(DP_2012);
indexer = DP2012(:,2);
fill = zeros(140,1);
select = rand(140,1)*300;
select = select + 1;
for i = 1:length(fill)
    fill(i) = round(select(i));    
end

for i = 1:140
    DP2012(9+i:150,2) = DP2012(fill(i),2);  
end

for i = 1:length(DP2012)
    if DP2012(i,2) <520
        DP2012(i,2) = 621;
    else
    end
end

figure;
x2 = DP2012(:,2);
y2 = DP2012(:,1);
plot(x2,y2)
%ylim([0 100])
xlabel('Density(kg m^{⁻3})');
ylabel('Depth(m)');
title('Random density profile')
set(gca, 'YDir','reverse')