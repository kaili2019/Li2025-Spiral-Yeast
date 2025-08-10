% this function get the absolute difference in area of a colony given it x 
% and y positions compared to target area.
%
% Kai Li
% 31 Mar 2025

function area = get_area(tt,xF,yF,target_area)
    x = xF(1:tt,:);
    y = yF(1:tt,:);
    
    % save colony as matrix
    [I_hist,~,~] = histcounts2(x+600, y+800,0:2:1200,0:2:1600);
    I = imresize(I_hist,2)>0.5;
    I = ~bwareaopen(~I, 10); % fill gaps in colony that has 10 or less pixels
    area = abs(sum(I(:))-target_area);
end