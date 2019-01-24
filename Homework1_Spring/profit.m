function pi= profit(s,svec,p)
% This function returns back the profit function of a firm
% If computed x is zero then it returns back -inf

   N = length(svec);
   Z = length(p);
   x = s*ones(N,1) - svec';
   
   pi = kron(p,x) - kron(ones(1,Z),0.2*x.^1.5);
   
   for i=1:N
       if x(i) < 0 
           pi(i,:) = -inf;
       end
   end
   


end

