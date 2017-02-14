function [Ix,Iy] = forward_gradient_2d(I)
%     Iy = volshift(I,1,0,0) - I;
%     Ix = volshift(I,0,1,0) - I;
%     Iz = volshift(I,0,0,1) - I;
%     Ix(:,n,:) = 0;
%     Iy(m,:,:) = 0;
%     Iz(:,:,o) = 0;
      [m,n] = size(I);
      Ix = zeros(size(I));
      Iy = Ix; 
          for y = 1:m-1
              for x = 1:n-1
                Ix(y,x) = I(y,x+1) - I(y,x);
                Iy(y,x) = I(y+1,x) - I(y,x);
              end
          end

      % if y = ymax;
  
          for x = 1:n-1
                Ix(m,x) = I(m,x+1) - I(m,x);
                Iy(m,x) = 0.0;   
          end
      
      %if x = xmax;
 
          for y = 1:m-1
                Ix(y,n) = 0.0;
                Iy(y,n) = I(y+1,n) - I(y,n);
          end
          
       Ix(m,n) = 0.0;
       Iy(m,n) = 0.0;
      
       Ix = single(Ix);
       Iy = single(Iy);
end