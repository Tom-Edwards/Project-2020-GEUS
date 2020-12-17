%% Indexes for each perturbation

clear
clc
thick_index = 10:2:30;
homo_index = 6:2:20;
homo_index = [homo_index 25 30 35];
ind1 = homo_index;
ind2 = thick_index;
fpath = 'toms';
PTT = 'Homo';

% Automate Allinone data processor

% read first index

for i = 1:length(ind1)

    
    PT1 = ind1(i);
run Allinone_new.m

end