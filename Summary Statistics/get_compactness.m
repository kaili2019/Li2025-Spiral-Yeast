% this function take a binary image I as input and computes it's
% compactness
% 
% Kai Li
% 21/10/2024

function compactness = get_compactness(I)
    
    area = sum(I(:));

    % get perimeter of image
    perimeterize = bwperim(I);

    perimeter = sum(perimeterize(:));

    compactness = (perimeter^2)/(area*4*pi);
end