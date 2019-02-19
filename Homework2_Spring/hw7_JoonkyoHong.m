% Code for HW2 in Spring (HW7) ECON512
% Written by Joonkyo Hong
% Build upon codes contained in code_QLD in ECON512 lecture
% Feb 16, 2019
% all is good, check plus
clear;

%% Question 1. Solve for the equilbrium through policy and value function iteration

run prmtr.m
global L eta kappa l nu delta betta lambda eps cost trans0 trans1
iter=0;
diff=10;
maxiter = 5000;

P0=repmat(cost,1,L);                      % competitive price
V0=zeros(L,L);


 tic
 while diff > eps && iter < maxiter
   
     % Policy function updating
           
           g = @(x) P_FOC(V0,x);
           P1 = fsolve(g,P0);
           
      % With updated P1, update the value function
      
           V1 = updateV(V0, P1, P0);
           
      % difference is
      
           diff = max(max(max(abs((V1-V0)./(1+V1)))),  max(max(abs((P1-P0)./(1+P1))))); 
           
      % Dampening
         
           V0 = (1-lambda)*V0 + lambda*V1;
           P0 = (1-lambda)*P0 + lambda*P1;
           
      % Go Next Iteration
      
       iter = iter+1;
       
       clc;
       disp('current iteration');
       disp(num2str(iter));
       disp('current difference');
       disp(num2str(diff));
       
       
 end   
 toc
 
 
 figure(1)
 surf((1:1:L),(1:1:L),V1);
 xlabel('w1');
 ylabel('w2');
 zlabel('value of the player 1');
 title("Value Function");
 
 figure(2)
 surf((1:1:L),(1:1:L),P1);
  xlabel('w1');
 ylabel('w2');
 zlabel('price charged by the player 1');
 title("Policy (Price) Function"); 
 
 
 %% Question 2. Compute the distribution of the state as time evolves
 
    % D0, D1, D2 in the equilibrium computed in Section 1.
    [D0eqm, D1eqm, D2eqm] = computeD(P1,P1');
    
    D0state = repmat(reshape(D0eqm',L*L,1),1,L*L);
    D1state = repmat(reshape(D1eqm',L*L,1),1,L*L);
    D2state = repmat(reshape(D2eqm',L*L,1),1,L*L); 
    
    % Conditional on each learning status, compute the transition matrix
    State0 = kron(trans0, trans0);
    State1 = kron(trans1, trans0);
    State2 = kron(trans0, trans1);
    
    % Finally the transition martrix is
    Pi = D0state.*State0 + D1state.*State1 + D2state.*State2;
    Pi = Pi./repmat(sum(Pi,2),1,L*L);     % To rule out the numerical noises making the sum exceed one
    
    % Compute the distribution of states after 10, 20, and 30 periods
    
    Pi10 = Pi^10;
    Pi20 = Pi^20;
    Pi30 = Pi^30;
    start = [1 zeros(1,899)];
   
    state10 = start*Pi10;
    state10 = reshape(state10,L,L);    % distribution of states, 10 periods later, L by L matrix
    state20 = start*Pi20;
    state20 = reshape(state20,L,L);    % distribution of states, 20 periods later, L by L matrix
    state30 = start*Pi30;
    state30 = reshape(state30,L,L);    % distribution of states, 30 periods later, L by L matrix
    
    figure(3)
    surf(1:1:L,1:1:L,state10);
    xlabel('w1');
    ylabel('w2');
    zlabel('probability mass');
    title('dist. of state, 10 periods later');
    
    figure(4)
    surf(1:1:L,1:1:L,state20);
    xlabel('w1');
    ylabel('w2');
    zlabel('probability mass');    
    title('dist. of state, 20 periods later');

    figure(5)
    surf(1:1:L,1:1:L,state30);
    xlabel('w1');
    ylabel('w2');
    zlabel('probability mass'); 
    title('dist. of state, 30 periods later');

    
%% Question 3. Compute the stationary distribution of the state

   start = [1 zeros(1,899)];

         diff=1;

         while diff>10e-16
               update=start*Pi;
               diff=norm(update-start)/norm(start);
               start=update;
         end
     
   stationary = reshape(update,L,L);
   
    figure(6)
    surf(1:1:L,1:1:L,stationary);
    xlabel('w1');
    ylabel('w2');
    zlabel('probability mass');
    title('stationary dist. of state');