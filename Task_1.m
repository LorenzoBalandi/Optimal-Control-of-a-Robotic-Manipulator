%% Final Project - 2 DOF Robotic Manipulator
% Optimal Control 2021
% Group 18: Balandi, Ghinelli, Prandin, January 2022
% Task 1

close all; clear; clc
tic % start counting time

% Run the code to generate the 'dynamics' function (see generatedynamics.m)
% generatedynamics

%% Parameters definition

max_iters = 2e2;
tol = 1e-6;

tf = 30; % seconds

% Parameters:
params.dyn.dt = 1e-3;
params.dyn.mm1 = 2;         % kg
params.dyn.mm2 = 2;         % kg
params.dyn.gg = 9.81;       % m/s^2
params.dyn.ll1 = 1;         % m
params.dyn.ll2 = 1;         % m
params.dyn.rr1 = 0.5;       % m
params.dyn.rr2 = 0.5;       % m
params.dyn.J_iner1 = 0.5;   % kg*m^2
params.dyn.J_iner2 = 0.5;   % kg*m^2

dt = params.dyn.dt;
mm1 = params.dyn.mm1;
mm2 = params.dyn.mm2;
gg = params.dyn.gg;
ll1 = params.dyn.ll1;
ll2 = params.dyn.ll2;
rr1 = params.dyn.rr1;
rr2 = params.dyn.rr2;
J_iner1 = params.dyn.J_iner1;
J_iner2 = params.dyn.J_iner2;

% We define the weights of the matrix Q
wq1 = 100;
wq2 = 1;
wq3 = 100;
wq4 = 1;
params.cost.QQ = [wq1, 0, 0, 0; ...
                  0, wq2, 0, 0; ...
                  0, 0, wq3, 0; ...
                  0, 0, 0, wq4];
params.cost.QQf = [wq1, 0, 0, 0; ...
                  0, wq2, 0, 0; ...
                  0, 0, wq3, 0; ...
                  0, 0, 0, wq4];           
% We define the weights of the matrix R
wr1 = 0.0007;
wr2 = 0.0007;
params.cost.RR = [wr1, 0;...
                  0, wr2];

TT = tf/params.dyn.dt;

state_dim = 4;
input_dim = 2;

% Flags for Armijo
ARMIJO_flag = 1;
gamma_fix = 0.1;

fprintf("Parameters defined\n")

%% Reference
% We define the reference angles and velocities (xx_ref) and the
% corresponding reference inputs (uu_ref)

% The robot starts from downward position (stable equilibrium)
ref_deg_q1_i = -90; % initial
ref_deg_q2_i = 0; % initial
ref_deg_q1_f = -10; % final
ref_deg_q2_f = 70; % final

xx_ref = zeros(state_dim, TT);
uu_ref = zeros(input_dim, TT);

xx_ref(1,1:TT/2) = deg2rad(ref_deg_q1_i);
xx_ref(1,TT/2:end) = deg2rad(ref_deg_q1_f);
xx_ref(3,1:TT/2) = deg2rad(ref_deg_q2_i);
xx_ref(3,TT/2:end) = deg2rad(ref_deg_q2_f);

% plots
figure(1);
subplot(1,2,1)
%plot(xx_ref(1,:),'LineWidth',2);
plot(rad2deg(xx_ref(1,:)),'LineWidth',2);
grid on
title('\theta_1 reference');
ylabel('\theta (deg)');
xlabel('t');
subplot(1,2,2)
%plot(xx_ref(3,:),'LineWidth',2);
plot(rad2deg(xx_ref(3,:)),'LineWidth',2);
grid on
title('\theta_2 reference');
ylabel('\theta (deg)');
xlabel('t');

% u_ref is such that it balances the g(q_ref) term. In this way, it keeps
% the manipulator in the desired equilibrium position, with q_dot and
% q_doubledot equal to zero
uu_ref(1,1:TT/2) = (mm1*rr1+mm2*ll1)*gg*cos(xx_ref(1,1))+mm2*gg*rr2*cos(xx_ref(1,1)+xx_ref(3,1));
uu_ref(1,TT/2:end) = (mm1*rr1+mm2*ll1)*gg*cos(xx_ref(1,TT))+mm2*gg*rr2*cos(xx_ref(1,TT)+xx_ref(3,TT));
uu_ref(2, 1:TT/2) = mm2*gg*rr2*cos(xx_ref(1,1)+xx_ref(3,1));
uu_ref(2, TT/2:end) = mm2*gg*rr2*cos(xx_ref(1,TT)+xx_ref(3,TT));

figure(2);
subplot(1,2,1)
plot(uu_ref(1,:),'LineWidth',2);
grid on
ylabel('u (Nm)');
xlabel('t');
title('u_1 reference');
subplot(1,2,2)
plot(uu_ref(2,:),'LineWidth',2);
grid on
ylabel('u (Nm)');
xlabel('t');
title('u_2 reference');

fprintf("Reference defined\n")

%% Trajectory definition
% So far, we have defined the reference (step).
% Now we define the optimal trajectory to pass from the initial equilibrium
% configuration (xx_ref(1),uu_ref(1)) to the final one (xx_ref(2),uu_ref(2))

xx = zeros(state_dim, TT, max_iters);
uu = zeros(input_dim, TT, max_iters);

% xx(1,:,1) = deg2rad(ref_deg_q1_i); % initialize x
% xx(3,:,1) = deg2rad(ref_deg_q2_i);
% uu(1,:,1) = uu_ref(1,1);
% uu(2,:,1) = uu_ref(2,1);

% Initialize x and u using a PD+gravity compensation control scheme
[xx(:,:,1),uu(:,:,1)] = algorithm_initialization(xx_ref, params, 1); % first iteration of x and u

JJ = zeros(max_iters,1);
descent = zeros(max_iters,1);

fprintf('-*-*-*-*-*-\n');

kk = 1; % iteration index

for tt=1:TT-1
  % In this for loop we build the FIRST iteration of the cost based on
  % state and input from the initial instant to the second last one
  [cost_dummy, ~] = stage_cost(xx(:,tt,kk), uu(:,tt,kk), xx_ref(:,tt), uu_ref(:,tt), params);
  JJ(kk) = JJ(kk) + cost_dummy;
end

% Here we define the last value of the cost, exploiting the state at the
% last instant (TT)
[cost_dummy, ~] = term_cost(xx(:,TT,kk), xx_ref(:,TT), params);
JJ(kk) = JJ(kk) + cost_dummy;

% MAIN LOOP (runs until max_iters or until tolerance is reached)
for kk=1:max_iters-1 
 
  KK = zeros(input_dim,state_dim, TT);
  sigma = zeros(input_dim, TT);
  
  pp = zeros(state_dim, TT);
  PP = zeros(state_dim,state_dim, TT);
  
  % Initialization of the terms p and P that have to be used to perform the
  % control:
  % p_T = q_T = gradient of the terminal cost
  % P_T = Q_T = hessian of the terminal cost
  % The gradient and the hessian are given as second and third outputs by 
  % the function term_cost
  [~, pp(:,TT), PP(:,:,TT)] = term_cost(xx(:,TT,kk), xx_ref(:,TT), params);
  
%   Backward iteration
  for tt = TT-1:-1:1
    
    [~, fx, fu, pfxx, pfuu, pfux] = dynamics(xx(:,tt,kk), uu(:,tt,kk),pp(:,tt+1));
%     pfxx = zeros(state_dim, state_dim); 
%     pfuu = zeros(input_dim, input_dim);
%     pfux = zeros(input_dim, state_dim);
    [~, lx, lu, lxx, luu, lux] = stage_cost(xx(:,tt,kk), uu(:,tt,kk),xx_ref(:,tt), uu_ref(:,tt), params);

    % Compute gain and ff descent direction
    
    KK(:,:,tt) = -(luu + fu*PP(:,:,tt+1)*fu'+ pfuu)\(lux + fu*PP(:,:,tt+1)*fx' + pfux);
    if (KK(:,:,tt)>1e5)
        fprintf("\nK > 1e5!!!\n")
        disp(tt)
        break;
    end
                              
    sigma(:,tt) = -(luu + fu*PP(:,:,tt+1)*fu'+ pfuu)\(lu + fu*pp(:,tt+1));
    if (sigma(:,tt)>1e5)
        fprintf("\nSigma > 1e5!!!\n")
        disp(tt)
        break;
    end                         
                              
    % Update PP and pp
    PP(:,:,tt) = (lxx + fx*PP(:, :, tt+1)*fx' + pfxx) - KK(:, :,tt)'*(luu + fu*PP(:,:,tt+1)*fu' + pfuu)*KK(:, :, tt);
    if (PP(:,:,tt)>1e5)
        fprintf("\nP > 1e5!!!\n")
        disp(tt)
        break;
    end 
  
    pp(:,tt) = (lx + fx*pp(:,tt+1))- KK(:, :,tt)'*(luu + fu*PP(:,:,tt+1)*fu' + pfuu)*sigma(:, tt);
    if (pp(:,tt)>1e5)
        fprintf("\np > 1e5!!!\n")
        disp(tt)
        break;
    end 

    descent(kk) = descent(kk) - sigma(:,tt)'*sigma(:,tt);
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %   ARMIJO gamma_stepsize selection
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %gamma_stepsize selection
  if ARMIJO_flag
    
    cc = 0.0001;
    rho = 0.5;
    gammas = 1;
    
    cost_arm = [];
    
    xx_temp = zeros(state_dim,TT);
    uu_temp = zeros(input_dim,TT);

    xx_temp(:,1) = xx_ref(:,1);
    JJtemp = 0;
    
    for tt = 1:TT-1
      
      
      uu_temp(:,tt) = uu(:,tt,kk) + gammas(end)*sigma(:,tt) + ...
        KK(:,:,tt)*(xx_temp(:,tt) - xx(:,tt,kk));
      
      [xx_temp(:,tt+1),~] = dynamics(xx_temp(:,tt),uu_temp(:,tt),...
        pp(:,tt+1));

      [cost_dummy, ~,~] = stage_cost(xx_temp(:,tt), uu_temp(:,tt), ...
        xx_ref(:,tt), uu_ref(:,tt), params);
      JJtemp = JJtemp + cost_dummy;
    end
    
    [cost_dummy, ~] = term_cost(xx_temp(:,TT), xx_ref(:,TT), params);
    
    JJtemp = JJtemp + cost_dummy;
    
    cost_arm = [cost_arm; JJtemp];
    
    % ARMIJO LOOP
    while cost_arm(end) > JJ(kk) + cc*gammas(end)*descent(kk)
      
      gammas = [gammas; gammas(end)*rho];
      
      % Evaluate cost for gamma_i
      xx_temp(:,1) = xx_ref(:,1);
      
      JJtemp = 0;
      
      for tt = 1:TT-1
        % Compute input
        uu_temp(:,tt) = uu(:,tt,kk) + gammas(end)*sigma(:,tt) + KK(:,:,tt)*(xx_temp(:,tt) - xx(:,tt,kk));
        %
        [xx_temp(:,tt+1),~] = dynamics(xx_temp(:,tt),uu_temp(:,tt), pp(:,tt+1));
        
        [cost_dummy, ~,~] = stage_cost(xx_temp(:,tt), uu_temp(:,tt), xx_ref(:,tt), uu_ref(:,tt), params);
        JJtemp = JJtemp + cost_dummy;
      end
      
      [cost_dummy, ~] = term_cost(xx_temp(:,TT), xx_ref(:,TT), params);
      JJtemp = JJtemp + cost_dummy;
      
      cost_arm = [cost_arm; JJtemp];
      
    end
    
    gamma_steps = gammas;
    gamma = gammas(end);
    
  else
    
    gamma = gamma_fix;
    
  end
  
  
 % Update trajectory
  
xx(:,1,kk+1) = xx_ref(:,1);
  for tt=1:TT-1

    uu(:,tt,kk+1) = uu(:,tt,kk) + gamma*sigma(:,tt) + KK(:,:,tt)*(xx(:,tt,kk+1) - xx(:,tt,kk));
    
    [xx(:,tt+1, kk+1),~] = dynamics(xx(:,tt,kk+1), uu(:,tt,kk+1), pp(:,tt+1));
 
    [cost_dummy, ~,] = stage_cost(xx(:,tt,kk+1), uu(:,tt,kk+1), xx_ref(:,tt), uu_ref(:,tt), params);
    
    JJ(kk+1) = JJ(kk+1) + cost_dummy;
    
  end

  [cost_dummy, ~] = term_cost(xx(:,TT,kk+1), xx_ref(:,TT), params);
  JJ(kk+1) = JJ(kk+1) + cost_dummy;
  
  fprintf('Iter: %d\n',kk);
  fprintf('descent: %.4e\n', descent(kk));
  fprintf('cost: %.4e\n', JJ(kk));
  
  if abs(descent(kk))<tol
    max_iters = kk;
    fprintf('Tolerance reached!\n');
    break;
  end
  
end 
% main loop

fprintf('The whole algorithm took %f seconds to run.\n',toc);

%% Add last samples (for plots)

uu(:,TT,max_iters) = uu(:,TT-1,max_iters);

%% Plots

% PLOT OF RESULTING STATE AND INPUT TRAJECTORIES

star = max_iters;

figure(3);
stairs(1:TT, xx(1,:,star),'LineWidth',2);
hold on;
stairs(1:TT, xx_ref(1,:),'--','LineWidth',2);
ylabel('x1_t (rad)');
xlabel('t');
grid on;
zoom on;

figure(4);
stairs(1:TT, xx(2,:,star),'LineWidth',2);
hold on;
stairs(1:TT, xx_ref(2,:),'--','LineWidth',2);
ylabel('x2_t (rad/s)');
xlabel('t');
grid on;
zoom on;

figure(5);
stairs(1:TT, xx(3,:,star),'LineWidth',2);
hold on;
stairs(1:TT, xx_ref(3,:),'--','LineWidth',2);
ylabel('x3_t (rad)');
xlabel('t');
grid on;
zoom on;

figure(6);
stairs(1:TT, xx(4,:,star),'LineWidth',2);
hold on;
stairs(1:TT, xx_ref(4,:),'--','LineWidth',2);
ylabel('x4_t (rad/s)');
xlabel('t');
grid on;
zoom on;

figure(7);
stairs(1:TT, uu(1,:,star),'LineWidth',2);
hold on;
stairs(1:TT, uu_ref(1,:),'--','LineWidth',2);
ylabel('u1_t (Nm)');
xlabel('t');
grid on;
zoom on;

figure(8);
stairs(1:TT, uu(2,:,star),'LineWidth',2);
hold on;
stairs(1:TT, uu_ref(2,:),'--','LineWidth',2);
ylabel('u2_t (Nm)');
xlabel('t');
grid on;
zoom on;

% Descent direction
figure(9);
semilogy(1:max_iters, abs(descent(1:max_iters)), 'LineWidth',2);
ylabel('descent');
xlabel('iter')
grid on;
zoom on;

% Cost error (normalized)
figure(10);
semilogy(1:max_iters, abs((JJ(1:max_iters)-JJ(max_iters))/JJ(max_iters)), 'LineWidth',2);
ylabel('J(u^k)-J(u^{max})/J(u^{max})');
xlabel('iter')
grid on;
zoom on;
