function z = indic_fcn(x,y)
% Indicator function 1(x^2 + y^2 <=1)

% x : Nx by 1 vector
% y : Ny by 1 vector
% This function will return z which is Nx by Ny matrix. 
% The (i,j)-th element of z will be
%               1,  if  x(i)^2 + y(j)^2 <= 1
%               0,  otherwise

 
  z = zeros(length(x),length(y));
  
  for i=1:length(x)
      for j=1:length(y)
      
            if x(i)^2+y(j)^2 > 1
                z(i,j) = 0;
            else
                z(i,j) = 1;
            end
      end
  end
       
end

