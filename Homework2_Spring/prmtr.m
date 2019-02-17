% SET UP THE PARAMETERS TO BE USED 
% Written by Joonkyo Hong 
% 16 FEB 2019

global L eta kappa l nu delta betta lambda eps cost trans0 trans1

L = 30;                      % know-how
eta = log(0.85)/log(2);      % slope of learning cost
kappa = 10;      % marginal cost of learning
delta = 0.03;    % depreciation of know-how
l     = 15;      % the units of know-how at that the learning curve flattens out
nu    = 10;      % quality of product
betta  = 1/(1.05); % discount factor

lambda = 1;     % dampening parameter
eps    = 10e-5;   % criterion


%%% Creating learning cost function
cost = zeros(L,1);
cost(1:l-1) = kappa*(1:1:l-1).^(eta);
cost(l:L)   = kappa*(l)^(eta);

%%% Transition Prob

   % when q=0
      trans0 = diag([1 (1-delta).^(2:1:L)],0) + diag(1-(1-delta).^(2:1:L),-1) ;
   % when q=1
      trans1 = diag([1-(1-delta).^(1:1:L-1) 1],0) + diag((1-delta).^(1:1:L-1),1);
      
