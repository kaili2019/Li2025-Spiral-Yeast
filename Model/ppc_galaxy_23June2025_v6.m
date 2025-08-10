% Implement off-lattice model using for data into Mixture Density Network.
% 
% 
% Kai Li
% 18 Dec 2024

clear 
clc
close all

rng(1102345) % set random seed for HPC use 

% Initialize theta_1
Nsims = 50; % number of iterations

samples = cell(Nsims,3);

% Total number of nutrients available 
N = 60000;

parameters = readmatrix("samples_23Jun2025_galaxy_r3_v2.csv");
parameters(:,6) = parameters(:,6)*pi/20;
idx = randi(length(parameters),Nsims,1);
sub_samples = parameters(idx,:);

%%
tStart = tic;
for i = 1:Nsims
    tic
    
    % MC3. Generate a data set x^* ~ f(x|thetaStar) 
    Telong = sub_samples(i,1);
    p2sProb = sub_samples(i,2);
    s2pProb = sub_samples(i,3); 
    pc = sub_samples(i,4); 
    pa = sub_samples(i,5); 
    angle_prolif = sub_samples(i,6);

    [file_name,x,y] = run_off_lattice_galaxy(N,Telong,p2sProb,s2pProb,pc,pa,angle_prolif);
        
    % save information from chain 
    samples{i,1} = x;
    samples{i,2} = y;
    samples{i,3} = [Telong,p2sProb,s2pProb,pc,pa,angle_prolif];

    fprintf("%d simulated completed and took %0.2f seconds with %d simulations left.\n", i-1, toc, Nsims-i);
    fprintf("Telong = %0.3f, p2sProb = %0.3f, s2pProb = %0.3f, pc = %0.3f, pa = %0.3f, angle = %0.3f\n\n",...
            Telong,p2sProb,s2pProb,pc,pa,angle_prolif);

    if mod(i,25) == 0 
        save("simulations/samples v1 i="+i+" "+datestr(datetime(now,'ConvertFrom','datenum'))+".mat","samples")
    end
    
end

tEnd = toc(tStart);

fprintf("total time taken is %0.2f seconds\n",tEnd);

%%
save("final_samples "+datestr(datetime(now,'ConvertFrom','datenum'))+".mat","samples")
