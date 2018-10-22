%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Homework #4 ECON 512                                    %
% Written by Joonkyo (Jay) Hong, 20 Oct 2018              %
% Modified on 22 Oct 2018                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;

addpath('./CEtools/');   % First, add path CEtools %
N = 100;                 % # of draws or nodes 

%% Questions 1 (Dart-throwing method with quasi-MC approach)

     [x1, ~]= qnwequi(N, [0 0], [1 1]);
     z   =  indic_fcn(x1(:,1),x1(:,2));
     
     pi1 =  4*mean(mean(z));
     error1 = abs(pi1 - pi);
     
%% Question 2 (Dart-throwing method with Newton-Coates approach)
 
     pi2 = 4*Int_indic([0 0],[1 1],N,N);     
     error2 = abs(pi2 - pi);

%% Questions 3 (Pythagorean method with quasi-MC approach)

     [x3, ~]  = qnwequi(N,0,1);
     pi3 = 4*mean(sqrt(1-x3.^2));
     error3 = abs(pi3 - pi);
     
%% Questions 4 (Pythagorean method with Newton-Coates approach)

     pi4 = 4*Int_simp(@(x) sqrt(1-x.^2), 0, 1, N);
     error4 = abs(pi4 - pi);
     
     
 result_from_1_to_4 =  ...
  [" ",           "approximate pi", "absolute error";
   "q-MC with DT",     pi1,               error1      ;
   "NC with DT"  ,     pi2,               error2      ;
   "q-MC with PYT",    pi3,               error3      ;
   "NC with PYT" ,     pi4,               error4      ];
%% Question 5 (Pseudo-MC)

N       = [100; 1000; 10000]; 
num_sim = 200;                   % # of simulations
store_array = zeros(2,length(N),num_sim);  % 3D-array that will store the Squared errors


for n=1:length(N)   
            
    for i=1:num_sim
        
          % Dart-Throwing with Pseudo-MC
        
          x1 = rand(N(n),2);    % Pseudo-MC 
         
          z   =  indic_fcn(x1(:,1),x1(:,2));
          pi1 =  4*mean(mean(z));
          store_array(1,n,i) = (pi1 - pi)^2;        
                    
          % Pythagorean with Pseudo-MC
        
          x3 = rand(N(n),1);    % Pseudo-MC 
         
          pi3 = 4*mean(sqrt(1-x3.^2));
          store_array(3,n,i) = (pi3 - pi)^2;
          
    end
    
          % Dart-Throwing with Newton-Cotes
           
          pi2 = 4*Int_indic([0 0],[1 1],N(n),N(n));     
          store_array(2,n,:) = (pi2 - pi)^2;
          
          % Pythagorean with Newton-Cotes

          pi4 = 4*Int_simp(@(x) sqrt(1-x.^2), 0, 1, N(n));
          store_array(4,n,:) = (pi4-pi)^2;
          
          
end



                   
results_mat = mean(store_array,3);   % Calculate mean over the simulations
results_mat= [[" ", "N=100", "N=1,000", "N=10,000"];
 ["pseudo-MC with DT"; "NC with DT"; "pseudo-MC with PYT"; "q-MC with PYT"], (results_mat)];

disp(" ");
disp(" ");
disp("Question 1~4. Approximated pi and abs error: N=100");
disp(result_from_1_to_4);

disp(" ");
disp("Question 5. Mean Squared Errors");
disp(results_mat);




      

