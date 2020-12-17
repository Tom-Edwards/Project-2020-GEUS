% Random spike generator - Create a random number of spike in a random
% number of places in the top 1000 positions (10m) of the site J profile

slab = ones(1000,1);                        % Empty slab
rl = round(100.*rand(1,1));                 % How many random lenses?
index = round((1000-10).*rand(rl,1) + 10);  % Positions for each lense
Low_limit = (500-250).*rand(1000,1) + 250;  % Lower limit for lense density
slab(index(:)) = Low_limit;                 % Paste into empty slab






% Now we need to create a density between the random lower limit and
% the upper limit in a for loop
lenses = zeros(1000,1);
for i = index(:):length(index)
    lenses(i) = (Low_limit(i)).*rand(1,1) + 350;      % Generate a random Density for lenses
end

         


%slab(index(:)) = lenses;                % index for existing profile
%plot(lenses)