function [v1,v2] = ev2d(grad,rho,Sx,Sy,um,vm,lambda,theta)
    [m,n] = size(um);
    v1 = zeros(m,n);
    v2 = v1;
    l_t =lambda*theta;
        for y = 1:m
            for x = 1:n
                dx = 0;dy = 0;
                rhoc = rho(y,x) + ( Sx(y,x).*um(y,x)+Sy(y,x).*vm(y,x));
                if rhoc < - l_t * grad(y,x)
                    dx = l_t * Sx(y,x);
                    dy = l_t * Sy(y,x);
                   
                else if rhoc > l_t * grad(y,x)
                        dx = -l_t * Sx(y,x);
                        dy = -l_t * Sy(y,x);
                       
                    else if grad(y,x) > 0
                          dx = - rhoc/grad(y,x) * Sx(y,x);
                          dy = - rhoc/grad(y,x) * Sy(y,x);
                         
                        end
                    end
                end
                v1(y,x) = um(y,x)+dx;
                v2(y,x) = vm(y,x)+dy;
            end
        end
end