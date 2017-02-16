function [D] = mydiv2d(PX,PY)
    [m,n] = size(PX);
    D = zeros(size(PX));
    L = zeros(size(PX));
   
        for yi = 2:m
            for xi = 2:n
               dvx = PX(yi,xi) - PX(yi,xi-1);
               dvy = PY(yi,xi) - PY(yi-1,xi);
         
               D(yi,xi) = dvx + dvy ;
               L(yi,xi) = L(yi,xi)+1;
            end
        end
    
    
    % y = 1;

        for xi = 2:n
            D(1,xi) = PX(1,xi) - PX(1,xi-1) + PY(1,xi);
            L(1,xi) = L(1,xi)+1;
        end

    
    % x = 1;
   
        for yi = 2:m
            D(yi,1) = PX(yi,1) + PY(yi,1) - PY(yi-1,1) ;
            L(yi,1) = L(yi,1)+1;
        end
   
    D(1,1) = PX(1,1) + PY(1,1);
    L(1,1) = L(1,1)+1;
end