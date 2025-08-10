% this function find the minimum radius of a given binary colony. 
%
% Kai Li 
% 26 June 2023 

function [x0, y0, rmin, rmax, rmean] = get_radii_rmean(I)
    
    % remove any particles that may cause errors in finding the center
    CC = bwconncomp(I);
    numPixels = cellfun(@numel,CC.PixelIdxList);
    [~,idx] = max(numPixels);
    J = zeros(size(I,1),size(I,2));
    J(CC.PixelIdxList{idx}) = 1;

    I = J;

    measurements = regionprops(I, 'Centroid');

    Image_centre = measurements.Centroid;

    dist = bwdist(I);
    
    % calculate maximum radius
    [y_max,x_max]=find(dist==0);
    
    xc_max = x_max-Image_centre(1);
    yc_max = y_max-Image_centre(2);
    
    E_dist_max = sqrt(xc_max.^2+yc_max.^2);

    rmax = max(E_dist_max);
    
    % calculate minimum radius
    [y_min,x_min]=find(dist~=0);
    xc_min = x_min-Image_centre(1);
    yc_min = y_min-Image_centre(2);

    E_dist_min = sqrt(xc_min.^2+yc_min.^2);
    rmin = min(E_dist_min);

    % calculate mean field density 

    x0 = Image_centre(1);
    y0 = Image_centre(2);
    
    pI= bwperim(I);
    
    [pY,pX] = find(pI);
    
    rDist = sqrt( (pX-x0).^2 + (pY-y0).^2 );

    rmean = mean(rDist);

end