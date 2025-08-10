% this function calculates the filament/colony area ratio of a binary yeast 
% colony image I. 
%
% AUTHOR: KAI LI 
% DATE: 22 Mar 2024
% 
% INPUT: 
%       I: a binary colony of the simulation or exeriment
%
% OUTPUT:
%       
%       filament_norm_ratio: normalised ratio between filamentous area with
%       respect to the whole colony area. 

function [outer_area, total_area] = get_colony_area_rmean(I) 
    
    [x0,y0,~,~] = get_radii(I);

    pI= bwperim(I);
    
    [pY,pX] = find(pI);
    
    rDist = sqrt( (pX-x0).^2 + (pY-y0).^2 );

    rmean = mean(rDist);

    % crop colony based on minimum radius 
    M = zeros(size(I,1),size(I,2));
    M(round(y0),round(x0)) = 1;
    R = bwdist(M);
    T = R >= rmean;    
    colony_outer = T.*I;
    
    % get filamentous to colony area ratio 
    outer_area = sum(colony_outer(:));
    total_area = sum(I(:));
    
end