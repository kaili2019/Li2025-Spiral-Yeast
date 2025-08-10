% this function calculates the number of subbranches of a colony based on
% the skeletonized image. 
% 
% AUTHOR: KAI LI 
% DATE: 22 Mar 2024
% 
% INPUT: 
%       I: a binary colony of the simulation or exeriment
%
% OUTPUT:
%       
%       norm_subbranch: normalised subbranch count from the skeletonization
%       method. 

function [subbranch, branch_length] = get_branch_length_rmean(I)

    % get csr radius

    [x0, y0, ~, ~, rmean] = get_radii_rmean(I);
    
    % get a matrix to crop out the centre of the colony based on csr radius
    M = zeros(size(I,1),size(I,2));
    M(round(y0),round(x0)) = 1;
    R = bwdist(M);
    C = R >= rmean;
    
    % crop a quarter of the colony
    T = zeros(size(I,1),size(I,2));
    T(1:round(y0),1:round(x0)) = 1;  

    % skeletonize
    skeleton = bwmorph(I.*T.*C, 'skel', Inf);

    % calculate sub branch number 
    It = logical(skeleton);
    [i,~] = find(bwmorph(It,'endpoints'));
    
    subbranch = length(i);

    branch_length = sum(It(:));
end