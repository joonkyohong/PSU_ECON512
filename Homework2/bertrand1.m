function fval = bertrand1(pa,pb,v)
% Simple Betrand competition model
% Multivariate Logit demand

% v: a vector of perceived values of a product
% p: a vector of prices

   p = [pa;pb];

% Demand vector
   D = exp(v-p)./(1+sum(exp(v-p)));
   
% Equations characterizing equilibrium

   fval = ones(length(p),1) - (eye(length(p))-diag(D))*p;
   fval = fval(1);
   
end

