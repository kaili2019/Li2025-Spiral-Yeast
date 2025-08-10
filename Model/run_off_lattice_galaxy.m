% this function take a set of parameter values and simulates the colony to
% the same area as the experiment upto 5% difference. 
% 
% AUTHOR: KAI LI
% 
% Date: 22 Mar 2024
% 
% INPUT: 
%       N: run colony upto this amount of cells
%       Telong: Proportion of cells to be sated before switching on 
%               pseudohyphal cells probabilities
%       p2sProb: Probability of pseudohyphal cells transitioning to 
%               sated cells
%       s2pProb: Probability of sated cells transitioning to pseudohyphal 
%               cells
%       pc: Probability of second daughter cell being sated given current 
%               cell is pseudohyphal
%       pa: Probability of choosing a sated cell to proliferate
%       upr_area: upper bound of area within 5% of experimental area
%       lwr_area: lower bound of area within 5% of experimental area
%       exp_area: experimental colony pixel area
%
% OUTPUT:
%       I: a binary colony of the simulation
%       file_name: the file name associated to the simualtion including the
%       information related to the particular simulation


function [file_name,x,y] = run_off_lattice_galaxy(N,Telong,p2sProb,s2pProb,pc,pa,angle_prolif)

    % run simualtion 
    [a2,b2,b_el2,pos] = simulate_colony_galaxy(N,Telong,p2sProb,s2pProb,pc,pa,angle_prolif);
    
    % simulates colony with area correction 
    
    % calculate (x,y) coordinates of cells within colony
    [x,y] = calc_colony_coordinates(a2,b2,b_el2,pos,N);
    
    % saving colony as tiff
    file_name = "Nfinal="+N+"_N="+N+"_Telong="+Telong+"_p2sProb="+p2sProb+"_s2pProb"...
                +s2pProb+"_pc="+pc+"_pa="+pa; 
    
end