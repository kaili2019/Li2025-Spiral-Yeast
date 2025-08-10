% Implement off-lattice model using for data into Neural Likelihood
% Estimator
% 
% 
% Kai Li
% 18 Dec 2024

clear 
clc
close all

rng(1102) % set random seed for HPC use 

% Initialize theta_1
Nsims = 150; % number of iterations

samples = cell(Nsims,3);

% Total number of nutrients available 
N = 60000;

%%
tStart = tic;
for i = 1:Nsims
    tic
    
    % propose parameter randomly
    Telong = rand;
    p2sProb = rand;
    s2pProb = rand; 
    pc = rand; 
    pa = rand; 
    angle_prolif = rand*(pi/20);

    [file_name,x,y] = run_off_lattice_galaxy(N,Telong,p2sProb,s2pProb,pc,pa,angle_prolif);
        
    % save information from chain 
    samples{i,1} = x;
    samples{i,2} = y;
    samples{i,3} = [Telong,p2sProb,s2pProb,pc,pa,angle_prolif];

    fprintf("%d simulated completed and took %0.2f seconds with %d simulations left.\n", i-1, toc, Nsims-i);
    fprintf("Telong = %0.3f, p2sProb = %0.3f, s2pProb = %0.3f, pc = %0.3f, pa = %0.3f, angle = %0.3f\n\n",...
            Telong,p2sProb,s2pProb,pc,pa,angle_prolif);

    if mod(i,50) == 0 
        save("simulations/samples v1 i="+i+" "+datestr(datetime(now,'ConvertFrom','datenum'))+".mat","samples")
    end
    
end

tEnd = toc(tStart);

fprintf("total time taken is %0.2f seconds\n",tEnd);

%% save data
save("final_samples "+datestr(datetime(now,'ConvertFrom','datenum'))+".mat","samples")
