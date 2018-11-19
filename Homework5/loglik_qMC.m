function neglik = loglik_qMC(para,dta,k)
% Return negative log likelihood function
% Built for the model described in HW5 ECON512
% Built upon the code for Roberts and Tybout (1997)
% Written by Joonkyo (Jay) Hong
% 18 Nov 2018

% INPUTs:
%     x, z:      regressors (vector)
%     y:         dependent variable (2-column matrix)
%     para:      vecter of parameters
%     k:         number of points for Gaussian quadrature (1-by-2 matrix)
% OUTPUTs:
%     neglik:    negative value of log-likelihood function

% Parameters
   gamma = para(1);
   beta0 = para(2);
   sigbeta = para(3);
        
% Data

   x = dta(:,2);
   z = dta(:,3);
   y = dta(:,1);
   nt = length(y);
        
% Compute likelihood using Halton Norm Shuffle

    Norm = haltonNormShuffle(k,1,1);
    beta = beta0 + sqrt(sigbeta)*Norm;
    
Y = repmat(y,1,k);
w = repmat(gamma*z,1,k)+kron(beta,x);
phi = exp(w)./(1+exp(w));
phi = reshape(prod(reshape(phi.^Y.*(1-phi).^(1-Y),20,[])),100,[]);
neglik = -sum(log(mean(phi,2)));

end

