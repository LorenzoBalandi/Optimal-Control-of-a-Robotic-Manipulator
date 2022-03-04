%% Final Project - 2 DOF Robotic Manipulator
% Optimal Control 2021
% Group 18: Balandi, Ghinelli, Prandin, January 2022
% Tasks 2 and 3

close all; clear; clc

tic % start counting time

%% Parameters definition

max_iters = 2e2;
tol = 1e-6;

tf = 30; % seconds

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
wq1 = 1000;
wq2 = 10;
wq3 = 1000;
wq4 = 10;
params.cost.QQ = [wq1, 0, 0, 0; ...
                  0, wq2, 0, 0; ...
                  0, 0, wq3, 0; ...
                  0, 0, 0, wq4];
params.cost.QQf = [wq1, 0, 0, 0; ...
                  0, wq2, 0, 0; ...
                  0, 0, wq3, 0; ...
                  0, 0, 0, wq4];         

% We define the weights of the matrix R
wr1 = 0.0006;
wr2 = 0.0006;
params.cost.RR = [wr1, 0;...
                  0, wr2];

TT = tf/params.dyn.dt;

state_dim = 4;
input_dim = 2;

% Flags for Armijo
visu_traj = 0;
plot_flag = 1;
ARMIJO_flag = 1;
gamma_fix = 0.5;

fprintf("Parameters defined\n")

%% Shape definition
% The workspace is the circle of radius 2 m around the origin

% Circle
r = 0.1; % m
%r = 0.2; % m

x_circle = zeros(24000,1); %round(2*pi/0.01)
y_circle = zeros(24000,1); %round(2*pi/0.01)
x_off = 1; % offet of the center of the circle
y_off = 1; % offet of the center of the circle

l1 = 1; % lengths of the joints
l2 = 1;

i = 1; % first index of a vector

for teta = 0:2.6180e-04:2*pi
    x_circle(i) = r*cos(teta)+x_off; % parametric expression of a circle
    y_circle(i) = r*sin(teta)+y_off;
    i = i+1;
end

figure(1)
plot(x_circle,y_circle,'LineWidth',2)
axis equal
grid on
hold on
plot(x_off,y_off,'.','LineWidth',2) % plot the cente of the circle
plot(0,0,'x','LineWidth',2) % plot the origin

%% Inverse kinematic
% Robotics book, p. 92

teta1 = zeros(size(x_circle)); % vectors of the same dimension of the number of points of the shaoe
teta2 = zeros(size(x_circle));

for i=1:size(x_circle,1)
    x = x_circle(i); % define the point for which we require the IK solution
    y = y_circle(i);
    % IK solution:
    cost2 = (x^2+y^2-l1^2-l2^2)/(2*l1*l2);
    sint2 = sqrt(1-cost2^2);

    if (imag(sint2)~=0) % check if point is outside workspace
        disp("IK ERROR");
    end

    teta2(i) = atan2(sint2,cost2);

    sint1 = ((l1+l2*cost2)*y-l2*sint2*x)/(x^2+y^2);
    cost1 = ((l1+l2*cost2)*x+l2*sint2*y)/(x^2+y^2);

    teta1(i) = atan2(sint1,cost1);
end

%% REFERENCE DEFINITION
% teta1 and teta2 are the new references to be followed with the DDP
% algorithm
% We must apply DDP to obtain the optimal trajectory, which is done
% minimizing the cost function determined by matrices Q, R and by the state
% and input references, which are teta1, teta2 and the torque balancing the
% gravity. The desired speed is, as always, zero 

% The reference is:
%   - teta = (-90,0) for 1 second;
%   - fifth order polynomial to the beginning of the circle for 5 seconds
%   - circle for the remaining 24 seconds

xx_ref = zeros(state_dim, TT);
uu_ref = zeros(input_dim, TT);

% Polynomial definition:
teta1_dot_init = (teta1(20)-teta1(1))/(20);
teta2_dot_init = (teta2(20)-teta2(1))/(20);

a01 = deg2rad(-90);
a11 = 0;
a21 = 0;
a31 = 1/(2*5001^3)*(20*(0.005013847243211-deg2rad(-90)) - 8*5001*teta1_dot_init);
a41 = 1/(2*5001^4)*(30*(-0.005013847243211+deg2rad(-90)) + 14*5001*teta1_dot_init);
a51 = 1/(2*5001^5)*(12*(0.005013847243211-deg2rad(-90)) - 6*5001*teta1_dot_init);
xx_ref(1,1:999) = deg2rad(-90);
for i=1000:6000
    % definition of the ramp
    % xx_ref(1,i) = ((0.005013847243211-deg2rad(-90))/5001)*(i-1000)+deg2rad(-90);

    % Definition of the polynomial trajectory
    xx_ref(1,i) = a01 + a31*(i-999)^3 + a41*(i-999)^4 + a51*(i-999)^5;
end
xx_ref(1,6001:end) = teta1(:);

a03 = 0;
a13 = 0;
a23 = 0;
a33 = 1/(2*5001^3)*(20*1.465602425754508 - 8*5001*teta2_dot_init);
a43 = 1/(2*5001^4)*(-30*1.465602425754508 + 14*5001*teta2_dot_init);
a53 = 1/(2*5001^5)*(12*1.465602425754508 - 6*5001*teta2_dot_init);
xx_ref(3,1:999) = deg2rad(0);
for i=1000:6000
    % definition of the ramp
    % xx_ref(3,i) = ((1.465602425754508-deg2rad(0))/5001)*(i-1000)+deg2rad(0);

    % Definition of the polynomial trajectory
    xx_ref(3,i) = a33*(i-999)^3 + a43*(i-999)^4 + a53*(i-999)^5;
end
xx_ref(3,6001:end) = teta2(:);


% u_ref is such that it balances the g(q_ref) term. In this way, it keeps
% the manipulator in the desired equilibrium position, with q_dot and
% q_doubledot equal to zero
for i=1:TT
    uu_ref(1,i) = (mm1*rr1+mm2*ll1)*gg*cos(xx_ref(1,i))+mm2*gg*rr2*cos(xx_ref(1,i)+xx_ref(3,i));
    uu_ref(2,i) = mm2*gg*rr2*cos(xx_ref(1,i)+xx_ref(3,i));
end

% plots
if plot_flag == 1
    figure(2);
    subplot(1,2,1)
    plot(rad2deg(xx_ref(1,:)),'LineWidth',2);
    grid on
    title('\theta_1 reference');
    ylabel('\theta (deg)');
    xlabel('t');
    subplot(1,2,2)
    plot(rad2deg(xx_ref(3,:)),'LineWidth',2);
    grid on
    title('\theta_2 reference');
    ylabel('\theta (deg)');
    xlabel('t');
    
    figure(3);
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
end

fprintf("Reference defined\n")

%% DDP IMPLEMENTATION
xx_star = zeros(state_dim, TT, max_iters);
uu_star = zeros(input_dim, TT, max_iters);

[xx_star(:,:,1), uu_star(:,:,1)] = algorithm_initialization(xx_ref, params, 2);

JJ = zeros(max_iters,1);
descent = zeros(max_iters,1);

fprintf('-*-*-*-*-*-\n');

kk = 1; % iteration index

for tt=1:TT-1
  % In this for loop we build the FIRST iteration of the cost based on
  % state and input from the initial instant to the second last one
  [cost_dummy, ~] = stage_cost(xx_star(:,tt,kk), uu_star(:,tt,kk), xx_ref(:,tt), uu_ref(:,tt), params);
  JJ(kk) = JJ(kk) + cost_dummy;
end

% Here we define the last value of the cost, exploiting the state at the
% last instant (TT)
[cost_dummy, ~] = term_cost(xx_star(:,TT,kk), xx_ref(:,TT), params);
JJ(kk) = JJ(kk) + cost_dummy;

% MAIN LOOP (runs until max_iters or until tolerance is reached)
for kk=1:max_iters-1 
  if visu_traj
    
    figure(30);
    stairs(1:TT, xx_star(1,:,kk),'LineWidth',2);
    hold on;
    stairs(1:TT, xx_ref(1,:),'--','LineWidth',2);
    ylabel('x1_t');
    xlabel('t');
    grid on;
    zoom on;
    hold off;
    drawnow

    figure(31);
    stairs(1:TT, xx_star(2,:,kk),'LineWidth',2);
    hold on;
    stairs(1:TT, xx_ref(2,:),'--','LineWidth',2);
    ylabel('x2_t');
    xlabel('t');
    grid on;
    zoom on;
    hold off;
    drawnow
 
    figure(32);
    stairs(1:TT, xx_star(3,:,kk),'LineWidth',2);
    hold on;
    stairs(1:TT, xx_ref(3,:),'--','LineWidth',2);
    ylabel('x3_t');
    xlabel('t');
    grid on;
    zoom on;
    hold off;
    drawnow
 
    figure(33);
    stairs(1:TT, xx_star(4,:,kk),'LineWidth',2);
    hold on;
    stairs(1:TT, xx_ref(4,:),'--','LineWidth',2);
    ylabel('x4_t');
    xlabel('t');
    grid on;
    zoom on;
    hold off;
    drawnow
    
    figure(34);
    stairs(1:TT, uu_star(1,:,kk),'LineWidth',2);
    hold on;
    stairs(1:TT, uu_ref(1,:),'--','LineWidth',2);
    ylabel('u1_t');
    xlabel('t');
    grid on;
    zoom on;
    hold off;
    drawnow

    figure(35);
    stairs(1:TT, uu_star(2,:,kk),'LineWidth',2);
    hold on;
    stairs(1:TT, uu_ref(2,:),'--','LineWidth',2);
    ylabel('u2_t');
    xlabel('t');
    grid on;
    zoom on;
    hold off;
    drawnow
   
    %pause % uncomment to stop at each iteration
    
  end
 
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
  [~, pp(:,TT), PP(:,:,TT)] = term_cost(xx_star(:,TT,kk), xx_ref(:,TT), params);
  
%   Backward iteration
  for tt = TT-1:-1:1
    
    [~, fx, fu, pfxx, pfuu, pfux] = dynamics(xx_star(:,tt,kk), uu_star(:,tt,kk),pp(:,tt+1));
    [~, lx, lu, lxx, luu, lux] = stage_cost(xx_star(:,tt,kk), uu_star(:,tt,kk),xx_ref(:,tt), uu_ref(:,tt), params);

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
      
      uu_temp(:,tt) = uu_star(:,tt,kk) + gammas(end)*sigma(:,tt) + ...
        KK(:,:,tt)*(xx_temp(:,tt) - xx_star(:,tt,kk));
      
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
        % Compute the input
        uu_temp(:,tt) = uu_star(:,tt,kk) + gammas(end)*sigma(:,tt) + KK(:,:,tt)*(xx_temp(:,tt) - xx_star(:,tt,kk));
        
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
  
xx_star(:,1,kk+1) = xx_ref(:,1);
  for tt=1:TT-1
    
    uu_star(:,tt,kk+1) = uu_star(:,tt,kk) + gamma*sigma(:,tt) + KK(:,:,tt)*(xx_star(:,tt,kk+1) - xx_star(:,tt,kk));
    
    [xx_star(:,tt+1, kk+1),~] = dynamics(xx_star(:,tt,kk+1), uu_star(:,tt,kk+1), pp(:,tt+1));

    [cost_dummy, ~,] = stage_cost(xx_star(:,tt,kk+1), uu_star(:,tt,kk+1), xx_ref(:,tt), uu_ref(:,tt), params);
    
    JJ(kk+1) = JJ(kk+1) + cost_dummy;
    
  end
  
  [cost_dummy, ~] = term_cost(xx_star(:,TT,kk+1), xx_ref(:,TT), params);
  JJ(kk+1) = JJ(kk+1) + cost_dummy;
  
  fprintf('Iter: %d\n',kk);
  fprintf('descent: %.4e\n', descent(kk));
  fprintf('cost: %.4e\n', JJ(kk));
  
  if abs(descent(kk))<tol
    max_iters = kk;
    fprintf('Tolerance reached!\n');
    break;
  end

end % main loop

%% Add last samples (for plots)

uu_star(:,TT,max_iters) = uu_star(:,TT-1,max_iters);

%% Plots

% PLOT OF RESULTING STATE AND INPUT TRAJECTORIES

star = max_iters;

if plot_flag == 1

    figure(4);
    stairs(1:TT, xx_star(1,:,star),'LineWidth',2);
    hold on;
    stairs(1:TT, xx_ref(1,:),'--','LineWidth',2);
    ylabel('x1_t (rad)');
    xlabel('t');
    grid on;
    zoom on;
    
    figure(5);
    stairs(1:TT, xx_star(2,:,star),'LineWidth',2);
    hold on;
    stairs(1:TT, xx_ref(2,:),'--','LineWidth',2);
    ylabel('x2_t (rad/s)');
    xlabel('t');
    grid on;
    zoom on;
    
    figure(6);
    stairs(1:TT, xx_star(3,:,star),'LineWidth',2);
    hold on;
    stairs(1:TT, xx_ref(3,:),'--','LineWidth',2);
    ylabel('x3_t (rad)');
    xlabel('t');
    grid on;
    zoom on;
    
    figure(7);
    stairs(1:TT, xx_star(4,:,star),'LineWidth',2);
    hold on;
    stairs(1:TT, xx_ref(4,:),'--','LineWidth',2);
    ylabel('x4_t (rad/s)');
    xlabel('t');
    grid on;
    zoom on;
    
    figure(8);
    stairs(1:TT, uu_star(1,:,star),'LineWidth',2);
    hold on;
    stairs(1:TT, uu_ref(1,:),'--','LineWidth',2);
    ylabel('u1_t (Nm)');
    xlabel('t');
    grid on;
    zoom on;
    
    figure(9);
    stairs(1:TT, uu_star(2,:,star),'LineWidth',2);
    hold on;
    stairs(1:TT, uu_ref(2,:),'--','LineWidth',2);
    ylabel('u2_t (Nm)');
    xlabel('t');
    grid on;
    zoom on;
    
    % Descent direction
    figure(10);
    semilogy(1:max_iters, abs(descent(1:max_iters)), 'LineWidth',2);
    ylabel('descent');
    xlabel('iter')
    grid on;
    zoom on;
    
    % Cost error (normalized)
    figure(11);
    semilogy(1:max_iters, abs((JJ(1:max_iters)-JJ(max_iters))/JJ(max_iters)), 'LineWidth',2);
    ylabel('J(u^k)-J(u^{max})/J(u^{max})');
    xlabel('iter')
    grid on;
    zoom on;

end

%% TASK 3
% In xx_star and uu_star are contained the optimal trajectories that our
% system has to follow in order to start from the vertical position and
% then draw the circle
% Now it's necessary to perform trajectory tracking, which is done
% solving a LQR optimal control problem. First of all, let's linearize the
% system about the optimal trajectory:
pp_useless = zeros(state_dim, TT);
x0 = [deg2rad(-90); 0; deg2rad(0); 0]; % Expected initialization
% x0 = [deg2rad(-23); 0; deg2rad(19); 0]; % Error in the initialization


% Linearization tensor. For every tt we have the linearization of the
% system in a certain instant.
AA_star = zeros(state_dim, state_dim, TT);
BB_star = zeros(state_dim, input_dim, TT);

for tt=1:TT
    [xxp_dummy, fx_star, fu_star, ~] = dynamics(xx_star(:,tt,star), uu_star(:,tt,star), pp_useless(:,tt));
    AA_star(:,:,tt) = fx_star';
    BB_star(:,:,tt) = fu_star';
end

% Now we have to define the weight matrices with wich we can perform the
% LQR control, considering that it's very important for us to precisely
% track the first and third component of the state

% If the initialization is the same of the previous task, the weights could
% be all ones. If the inizialization is different from the previous task or
% the dynamics parameters changes, these weights allow the robot to
% precisely track the trajectory.
wq_star1 = 1000;
wq_star2 = 10;
wq_star3 = 1000;
wq_star4 = 10;
QQ_star = [wq_star1, 0, 0, 0; ...
           0, wq_star2, 0, 0; ...
           0, 0, wq_star3, 0; ...
           0, 0, 0, wq_star4];
wr_star1 = 1;
wr_star2 = 1;
RR_star = [wr_star1, 0;...
           0, wr_star2];

% QQf_star = idare(AA_star(:,:,TT),BB_star(:,:,TT),QQ_star,RR_star,[],[]);
QQf_star = QQ_star;
xx = zeros(state_dim,TT);
xx(:,1) = x0;
uu = zeros(input_dim,TT);
PP = zeros(state_dim,state_dim,TT);
KK_star = zeros(input_dim,state_dim,TT);
PP(:,:,end) = QQf_star; 

% Backward iteration to compute P
for tt = TT-1:-1:1
    AAt = AA_star(:,:,tt);
    BBt = BB_star(:,:,tt);
    PPtp = PP(:,:,tt+1);        
    PP(:,:,tt) = QQ_star + AAt'*PPtp*AAt - (AAt'*PPtp*BBt)*inv(RR_star + BBt'*PPtp*BBt)*(BBt'*PPtp*AAt);
end

% Forward iteration to compute K,u,x_t+1
for tt = 1:TT-1
    AAt = AA_star(:,:,tt);
    BBt = BB_star(:,:,tt);
    PPtp = PP(:,:,tt+1);
    KK_star(:,:,tt) = -(RR_star + BBt'*PPtp*BBt)\(BBt'*PPtp*AAt);

    uu(:,tt) = uu_star(:,tt, star) + KK_star(:,:,tt)*(xx(:,tt)-xx_star(:,tt, star));
    [xx(:,tt+1), ~] = dynamics(xx(:,tt), uu(:,tt), pp_useless(:,tt));
end 

uu(:,TT) = uu(:,TT-1);

fprintf('The whole algorithm took %f seconds to run.\n',toc);

% Code for passing data to Simscape Multibody
u1 = zeros(TT,2);
u2 = zeros(TT,2);
x1 = zeros(TT,2);
x2 = zeros(TT,2);
x3 = zeros(TT,2);
x4 = zeros(TT,2);
for jj=0:TT-1
    x1(jj+1,1) = dt*jj;
    x1(jj+1,2) = xx(1,jj+1);
    x2(jj+1,1) = dt*jj;
    x2(jj+1,2) = xx(2,jj+1);
    x3(jj+1,1) = dt*jj;
    x3(jj+1,2) = xx(3,jj+1);
    x4(jj+1,1) = dt*jj;
    x4(jj+1,2) = xx(4,jj+1);
    u1(jj+1,1) = dt*jj;
    u1(jj+1,2) = uu(1,jj+1);
    u2(jj+1,1) = dt*jj;
    u2(jj+1,2) = uu(2,jj+1);
end

% From now on, there are the plots
figure(12);
stairs(1:TT, xx(1,:),'LineWidth',2);
hold on;
stairs(1:TT, xx_star(1,:,star),'--','LineWidth',2);
ylabel('x1_t (rad)');
xlabel('t');
grid on;
zoom on;

figure(13);
stairs(1:TT, xx(2,:),'LineWidth',2);
hold on;
stairs(1:TT, xx_star(2,:,star),'--','LineWidth',2);
ylabel('x2_t (rad/s)');
xlabel('t');
grid on;
zoom on;

figure(14);
stairs(1:TT, xx(3,:),'LineWidth',2);
hold on;
stairs(1:TT, xx_star(3,:,star),'--','LineWidth',2);
ylabel('x3_t (rad)');
xlabel('t');
grid on;
zoom on;

figure(15);
stairs(1:TT, xx(4,:),'LineWidth',2);
hold on;
stairs(1:TT, xx_star(4,:,star),'--','LineWidth',2);
ylabel('x4_t (rad/s)');
xlabel('t');
grid on;
zoom on;

figure(16);
stairs(1:TT, uu(1,:),'LineWidth',2);
hold on;
stairs(1:TT, uu_star(1,:,star),'--','LineWidth',2);
ylabel('u1_t (Nm)');
xlabel('t');
grid on;
zoom on;

figure(17);
stairs(1:TT, uu(2,:),'LineWidth',2);
hold on;
stairs(1:TT, uu_star(2,:,star),'--','LineWidth',2);
ylabel('u2_t (Nm)');
xlabel('t');
grid on;
zoom on;

%% Visualisation of the end-effector trajectory
% Run this section to visualize the end-effector trajectory from the
% Simscape simulation. Before running this section, it is necessary to run
% the Simscape Multibody model, otherwise no x_motion and y_motion exist.

% figure(17); 
% x_motion = out.x_motion;
% y_motion = out.y_motion;
% comet(x_motion, y_motion);
% %plot(x_motion, y_motion,'LineWidth',2);
% axis equal;
% xlabel('X Coordinate (m)'); 
% ylabel('Y Coordinate (m)'); 
% grid on;
