%% Depth convertor for rescaling the depth input to the output 
% ind makes an array showing all values between two depths e.g. [1 1 1 2 2 2 2 3 4 4 4]
% good_ind selects on inputs that are numbers
% Find index finds all NAN's, then finds the index of the first NAN value, -1 to
% take the preceding value before the last. 
% Then you have the length of the vector containing all values below 

% choose resolution of depth averager (1=per metre, 0.5=half metre etc)


function  rho_out = ConvertDepth(xi_interpolated_depths,yi_interpolated_densities,res)
    ind = discretize(xi_interpolated_depths,res);
    good_ind = 1:(find(isnan(ind),1,'first')-1);
    rho_out = accumarray(ind(good_ind), yi_interpolated_densities(good_ind)',[],@mean);
     
end


% figure
% stairs(depth_in1,rho_in1)
% hold on
% stairs(depth_out(1:end-1),rho_out)