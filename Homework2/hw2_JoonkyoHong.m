%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Homework #2 ECON 512                                    %
% Written by Joonkyo (Jay) Hong, 20 Sept 2018             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;

%% Question 1

v=[2;2];
p=[1;1];

demand = exp(v-p)./(1+sum(exp(v-p)));

disp('D_A         D_B')
disp(num2str(demand'));
disp("Question 1 is done. Press any key to continue");
pause
%% Question 2

f = @(p) bertrand(p,v); %make function handle whose argument is a price vector
                        %Bertrand is an outside function that characterize
                        %equilibrium conditions

% START BROYDEN METHODS

maxiter = 100;
tol = 10e-8;

% initial guess
p = [1;1];

init_f = f(p);
invJ = eye(length(p));

tic
 for i=1:maxiter
      fnorm = norm(init_f);
      fprintf('iter %d: p(1) = %f, p(2) = %f, norm(f(x)) = %.8f\n', i, p(1), p(2), fnorm);
      if fnorm < tol
          break;
      end
      
      d = - (invJ*init_f);
      p = p + d;
      new_f = f(p);
      u = invJ*(new_f - init_f);
      invJ = invJ + ( (d - u) * (d'*invJ) )/ (d'*u);
      
      init_f = new_f;
 end
t=toc;

disp("Equilibrium Prices")
disp('p_A         p_B')
disp(num2str(p'));
disp("Computation Time")
disp(num2str(t));
disp("Question 2 is done. Press any key to continue");
pause

%% Question 3

pa = 1; pb=1;

f = @(p) bertrand(p,v); % Again, make function handle
fval = f([pa;pb]);
pa_old = 0;
pb_old = 0;
bigloop_iter = 1000;

tic

for i=1:bigloop_iter % Big loop
    if norm(fval) < 10e-8    
        break
    end

% Sub loop for good A 
% In order to perform secent method, fix two initial conditions

      fa = @(pa) bertrand1(pa,pb,v);   % FOC only for good A
      faold = fa(pa_old);

      maxiter = 100;
      tol = 10e-8;

      for j=1:maxiter
          faval = fa(pa);
          if norm(faval) < tol
             break
          else
             pa_new = pa - ( (pa - pa_old)/(faval-faold) )*faval;
             pa_old = pa;
             pa = pa_new;
             faold = faval;
          end
      end

% Sub loop for good B

       fb = @(pb) bertrand2(pa,pb,v);
       fbold = fb(pb_old);

     for j=1:maxiter
          fbval = fb(pb);
         if norm(fbval) < tol
            break
         else
            pb_new = pb - ( (pb - pb_old)/(fbval-fbold) )*fbval;
            pb_old = pb;
            pb = pb_new;
            fbold = fbval;
         end
     end
 
     fval=f([pa;pb]);     % Updating fval for the next loop
     fprintf('iter %d: p(1) = %f, p(2) = %f, norm(f(x)) = %.8f\n',i, pa, pb, norm(fval));
     
end
t=toc;

disp("Equilibrium Prices")
disp('p_A         p_B')
disp(num2str([pa pb]));
disp("Computation Time")
disp(num2str(t));
disp("Question 3 is done. Press any key to continue");
pause
%% Question 4

% Updating Rule

g = @(p) 1./(1 - exp(v-p)./(1+sum(exp(v-p))));

p = [1;1];

maxit = 1000;

tic
 for i=1:maxit
     nextp = g(p);
     diff = norm(nextp-p);
     fprintf('iter %d: p(1) = %f, p(2) = %f, norm(f(x)) = %.8f\n',i, p(1), p(2), diff);
      if diff < 10e-8
          break;
      end
      p = nextp;
 end
t=toc;

% What if v2=3?

v = [2;3];
g = @(p) 1./(1 - exp(v-p)./(1+sum(exp(v-p))));

p2 = [1;2];

maxit = 1000;

tic
 for i=1:maxit
     nextp2 = g(p2);
     diff = norm(nextp2-p2);
     fprintf('iter %d: p(1) = %f, p(2) = %f, norm(f(x)) = %.8f\n',i, p2(1), p2(2), diff);
      if diff < 10e-8
          break;
      end
      p2 = nextp2;
 end
t2=toc;

disp("CASE 1: v1=v2=2")
disp("Equilibrium Prices")
disp('p_A         p_B')
disp(num2str(p'));
disp("Computation Time")
disp(num2str(t));
disp("Question 4 is done. Press any key to continue");

disp(" ")

disp("CASE 2: v1=2, v2=3")
disp("Equilibrium Prices")
disp('p_A         p_B')
disp(num2str(p2'));
disp("Computation Time")
disp(num2str(t2));
disp("Question 4 is done. Press any key to continue");

pause
%% Question 5

vb_vec = 0:0.2:3;                        % vector of values 
p_eqmvec = zeros(length(vb_vec),2);      % vector that will contain the equilibrium prices
p = [1;1];                               % Initial guess of solution 


for iter=1:length(vb_vec)   
    
    vb = vb_vec(iter);
    v = [2;vb];
    f = @(p) bertrand(p,v);    
    
          % Perform Broyden Method to solve the equilibrium
          
          init_f = f(p);
          invJ = eye(length(p)); 
          maxiter = 1000;
          tol = 10e-8;
              
             for i=1:maxiter
                 fnorm = norm(init_f);
                 
                 if fnorm < tol
                     break;
                 end
                 
                 d = - (invJ*init_f);
                 p = p + d;
                 new_f = f(p);
                 u = invJ*(new_f - init_f);
                 invJ = invJ + ( (d - u) * (d'*invJ) )/ (d'*u);
                 
                 init_f = new_f;
             end
            
        
        p_eqmvec(iter,:) = p';
    
end

figure(1)
plot(vb_vec,p_eqmvec);
xlabel('v_{B}');
ylabel('Prices');
legend('p_{A}', 'p_{B}');
disp("Question 5 is done. See the figure that poped up");

    