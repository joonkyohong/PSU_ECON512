%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Homework #5 ECON 512                                    %
% Written by Joonkyo (Jay) Hong, 18 Nov 2018              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;

addpath('./CEtools/');   % First, add path CEtools %

% Load dataset

load 'hw5.mat';

x=reshape(data.X,100*20,1);
z=reshape(data.Z,100*20,1);
y=reshape(data.Y,100*20,1);
dta = [y x z];

% Option set for fmincon
options_opt = optimoptions('fmincon', 'Algorithm', 'interior-point', 'Display', 'iter', 'StepTolerance',1e-30,...
      'TolFun',1e-4,'MaxIter',1e2, 'MaxFunEvals', 1e8);
  
% Initial Values for Q1 and Q2 (Logit Estimates)
init_theta = glmfit([z x],y,'binomial','link','logit','constant','off');
init_theta = [init_theta' 1];

% Initial Values for Q4 (Logit Estimates)
init_theta2 = glmfit([z x],y,'binomial','link','logit');
init_theta2 = [init_theta2(2:3,1)' init_theta2(1,1) 1 1 0];


%% Question 1. Simple with Gaussian Quadrature Approach

    k = 20;                                     % # of nodes for GQ approach
    lik_fcn1 = @(theta) loglik_GQ(theta,dta,k);

  % Constraints in parameter space
  % gamma, beta0, sigma_beta

   lb = [-inf -inf   0];    
   ub = [inf   inf inf];

  % Initial Values for Optimization and Initial value of likelihood function
  
   fcnval_for_Q1 = -lik_fcn1([0, 0.1, 1]);
   init_fcn1 = -lik_fcn1(init_theta);
   
  % Minimizing minus log likelihood

   tic
   [theta_mle1,fval1,exitflag1,~] = fmincon(lik_fcn1,init_theta,[],[],[],[],lb,ub,[],options_opt);
   time1 = toc;

%% Question 2. Simple with Quasi-Monte Carlo

    k = 100;                                     % # of nodes for q-MC approach
    lik_fcn2 = @(theta) loglik_qMC(theta,dta,k);

  % Constraints in parameter space
  % gamma, beta0, sigma_beta

   lb = [-inf -inf   0];    
   ub = [inf   inf inf];

  % Initial Values for Optimization and Initial value of likelihood function
  
   fcnval_for_Q2 = -lik_fcn2([0, 0.1, 1]);   
   init_fcn2 = -lik_fcn2(init_theta);
   
  % Minimizing minus log likelihood

   tic
   [theta_mle2,fval2,exitflag2,~] = fmincon(lik_fcn2,init_theta,[],[],[],[],lb,ub,[],options_opt);
   time2 = toc;


%% Question 4. Full with Quasi-Monte Carlo

    k = 200;                                     % # of nodes for q-MC approach
    lik_fcn4 = @(theta) loglik_full(theta,dta,k);

  % Constraints in parameter space
  % gamma, beta0, u0, sigma_beta, sigma_u, rho

   lb = [-inf -inf -inf   0   0 -1];    
   ub = [inf   inf  inf inf inf  1];

  % Initial Values for Optimization and Initial value of likelihood function
  
   init_fcn4 = -lik_fcn4(init_theta2);
   
  % Minimizing minus log likelihood

   tic
   [theta_mle4,fval4,exitflag4,~] = fmincon(lik_fcn4,init_theta2,[],[],[],[],lb,ub,[],options_opt);
   time4 = toc;


%% Question 5. Report the results


 disp("              ");
 disp("Case 1: Simple with Gaussian Quadrature");
 disp("Initial");
 disp(num2str(init_theta));
 disp("MLE Estimates");
 disp(num2str(theta_mle1));
 disp("Function Value");
 disp("Initial");
 disp(num2str(init_fcn1));
 disp("Maximand");
 disp(num2str(-1*fval1));


 disp("              ");
 disp("Case 2: Simple with quasi-Monte Carlo");
 disp("Initial");
 disp(num2str(init_theta));
 disp("MLE Estimates");
 disp(num2str(theta_mle2));
 disp("Function Value");
 disp("Initial");
 disp(num2str(init_fcn2));
 disp("Maximand");
 disp(num2str(-1*fval2));

 
 disp("              ");
 disp("Case 3: Full with quasi-Monte Carlo");
 disp("Initial");
 disp(num2str(init_theta2));
 disp("MLE Estimates");
 disp(num2str(theta_mle4));
 disp("Function Value");
 disp("Initial");
 disp(num2str(init_fcn4));
 disp("Maximand");
 disp(num2str(-1*fval4));
 