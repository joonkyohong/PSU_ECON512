% Code for HW1 ECON512
% Written by Joonkyo Hong
% Build upon stochgrow.m in ECON512 lecture
% Jan 23, 2019
% All is good, check plus

%% Section 1. Solve for Value functions

% Parameter values

 delta=0.95; % discount factor
  
 Z=21;          % grid number for shocks
 sigma=0.1;     % std of error term in AR(1) process
 mu=0.5;        % unconditional mean 
 ro=0.5;        % AR(1) parameter 
 
 

% GRID FOR STOCK

N=200;

s0 = 100;              % Initial value of stock 
slow = 0;    % In the end, the stock will be zero because there's no growth
grids = ((s0)-slow)/(N-1);
s = slow:(grids):(s0);

% GRID FOR PRICE PROCESS

% grid for shock by using Tauchen's method for finite state Markov
% approximation

[prob,grid]=tauchen(Z,mu,ro,sigma);
disp(['The dimensions of prob are ' num2str(size(prob)) ])
disp(['The dimensions of grid are ' num2str(size(grid)) ])

% VALUE FUNCTION ITERATION

vinitial=zeros(N,Z);
vrevised=zeros(N,Z);
decision=zeros(N,Z);
  
invest=kron(ones(1,Z),s');
disp(['The dimensions of invest  are ' num2str(size(invest)) ])
ONE=ones(N,1);
 
%iteration

diff=1;

while diff>0.001
    
    Ev=vinitial*prob';   % find the expected value of value function
      
    for i=1:N           % for each k, find the optimal decision rule
      pi = profit(s(i),s,grid);  
      [vrevised(i,:),decision(i,:)]=max(pi+delta*Ev);
    end
    
     diff=norm(vrevised-vinitial)/norm(vrevised)
     vinitial=vrevised;
    
end

% compute decision rule

derule=zeros(N,Z);

for j=1:Z
    derule(:,j)=s(decision(:,j))';
end

net_invest = kron(ones(1,Z),s') - derule;

figure(1)
plot(s,[vrevised(:,8) vrevised(:,11) vrevised(:,14)]);
xlabel('initial stock');
ylabel('values');
title('value function');


figure(2)
plot(s, [net_invest(:,8) net_invest(:,11) net_invest(:,14)]);
xlabel('current stock');
ylabel('harvest');
title('optimal harvest policy function');


%% Section 2. Predicting remaining amounts of stock

% Compute the transition probability of state

     P=zeros(Z*N,Z*N);

     for i=1:Z
         for j=1:Z
              P((i-1)*N+1:i*N,(j-1)*N+1:j*N)=prob(i,j)*(kron(ones(1,N),derule(:,i))==kron(ones(N,1),s));
        
         end
     end
     

% Initial probability for the state with stock 100 and price 1

    index_1 = find(grid==1);
    pinitial = zeros(1,N*Z);
    pinitial(1,index_1*N)=1;
    
    
% Simulate the model

    tau = 20;           % # of Periods for expectation
    predicted = zeros(tau,3);
    
    
   for t=1:20
       
       px = pinitial*P;
       probs=zeros(N,Z);
       
       for i=1:Z
           probs(:,i) = px((i-1)*N+1:i*N)';
       end
       
       probs = sum(probs');
       cdf = cumsum(probs);
       lower = max(find(cdf<=0.05));
       upper = min(find(cdf>=0.95));
       
       mean = probs*s';
       lowerbdd = s(lower);
       upperbdd = s(upper);
       
      predicted(t,:) = [lowerbdd mean upperbdd];
      
    
   pinitial = px;
      
   end
   
       
figure(3)
plot(0:1:20, [100 100 100;predicted]);
xlabel('Periods');
ylabel('Predicted Remainig Amounts of Stock');
title('Predicted trajectory of amounts of stock');