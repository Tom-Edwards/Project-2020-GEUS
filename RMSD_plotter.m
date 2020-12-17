%% Root Mean square plotter

% Can be formatted into cell array?

thick_index = 10:2:30;
homo_index = 6:2:20;

ind1 = homo_index;
ind2 = thick_index;

number = '10';
outputFP = append('/PHomo',number);
RMSD_path = append('./toms/PerbsWD89_19',outputFP,'/RMSD_2012_Finsim.txt');
RMSD_2012_Sim_10 = readtable(RMSD_path);

for i = ind1(1):ind1(end)

    N = num2str(i);
RMSD_2012 = append('RMSD_2012_Sim_',N);
RMSD_array = table2array(RMSD_2012);

end

