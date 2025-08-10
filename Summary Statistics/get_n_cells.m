% This function given a lower bound a and upper bound b, colony coordinates
% xF and yF, target area, it estimates the closest number of cells that
% will give the target area.
% This uses a Golden section search approach:
% https://en.wikipedia.org/wiki/Golden-section_search
%
% Kai Li
% 31 Mar 2025

function Ncells = get_n_cells(a,b,xF,yF,target_area)

    invphi = (sqrt(5)-1)/2; % 1/phi
    
    while b - a > 100
        
        c = b - (b-a) * invphi;
        d = a + (b-a) * invphi;
    
        if get_area(round(c),xF,yF,target_area) < get_area(round(d),xF,yF,target_area)
            b = d;
        else
            a = c;
        end
        
    end

    Ncells = round((a+b)/2);
end