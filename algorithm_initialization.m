function [output1,output2] = algorithm_initialization(xx_ref, params, task)
% Algorithm initialization
% Script for the definition of an initialization "near" the optimal one
% Industrial robotics control scheme (PD+gravity compensation)
% clear all; clc;

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

    TT = 30/dt;
    state_dim = 4;
    input_dim = 2;
    max_iters = 2e2;
    plot_flag = 1;
    
    % Reference:
%     ref_deg_q1_i = -90; % initial
%     ref_deg_q2_i = 0; % initial
%     ref_deg_q1_f = -10; % final
%     ref_deg_q2_f = 70; % final
    
%     xx_ref = zeros(state_dim, TT);
%     uu_ref = zeros(input_dim, TT);
%     xx_ref(1,1:TT/2) = deg2rad(ref_deg_q1_i);
%     xx_ref(1,TT/2:end) = deg2rad(ref_deg_q1_f);
%     xx_ref(3,1:TT/2) = deg2rad(ref_deg_q2_i);
%     xx_ref(3,TT/2:end) = deg2rad(ref_deg_q2_f);

    x_tilda = zeros(input_dim,TT); % state error

% Different control parameters depending on task 1 or 2
    if task == 1
        Kp = [5, 0;...
              0, 5];
        Kd = [10, 0;...
              0, 10];
    end

    if task == 2
        Kp = [45, 0;...
              0, 45];
        Kd = [25, 0;...
              0, 25];
    end

    gravity = zeros(input_dim,TT);
    kk = 1; % first iteration
    xx_init = zeros(state_dim, TT, max_iters);
    uu = zeros(input_dim, TT, max_iters);
    uu_init = zeros(input_dim, TT, max_iters);

    xx_init(1,:,1) = deg2rad(-90); % initialize x
    xx_init(3,:,1) = deg2rad(0);

    pp = zeros(state_dim,TT); % initialized to 0 because it is not needed
    xx_dot = zeros(state_dim/2,TT);

    % PD + gravity compensation control scheme
    for i=1:TT-1
        % Compute errors
        x_tilda(1,i) = xx_ref(1,i) - xx_init(1,i,kk);
        x_tilda(2,i) = xx_ref(3,i) - xx_init(3,i,kk);
        gravity(1,i) = (mm1*rr1+mm2*ll1)*gg*cos(xx_init(1,i,kk))+mm2*gg*rr2*cos(xx_init(1,i,kk)+xx_init(3,i,kk));
        gravity(2,i) = mm2*gg*rr2*cos(xx_init(1,i,kk)+xx_init(3,i,kk));
        xx_dot(1,i) = xx_init(2,i,1); 
        xx_dot(2,i) = xx_init(4,i,1);
        % Update u
        uu(:,i,kk) = gravity(:,i) + Kp*x_tilda(:,i) - Kd*xx_dot(:,i);
        % Update dynamics
        [xx_init(:,i+1, kk),~] = dynamics(xx_init(:,i,kk), uu(:,i,kk),pp(:,i));
    end
    % Update last samples
    gravity(1,TT) = (mm1*rr1+mm2*ll1)*gg*cos(xx_init(1,TT,kk))+mm2*gg*rr2*cos(xx_init(1,TT,kk)+xx_init(3,TT,kk));
    gravity(2,TT) = mm2*gg*rr2*cos(xx_init(1,TT,kk)+xx_init(3,TT,kk));
    uu_init(1,:,kk) = gravity(1,:);
    uu_init(2,:,kk) = gravity(2,:);
    output1 = xx_init(:,:,1);
    output2 = uu_init(:,:,1);

    if plot_flag    % plots
        
        figure(21);
        subplot(1,2,1)
        plot(rad2deg((xx_init(1,:,1))),'LineWidth',2);
        grid on
        title('\theta_1 init');
        ylabel('\theta (deg)');
        xlabel('t');
        subplot(1,2,2)
        plot(rad2deg((xx_init(3,:,1))),'LineWidth',2);
        grid on
        title('\theta_2 init');
        ylabel('\theta (deg)');
        xlabel('t');

        figure(22);
        subplot(1,2,1)
        plot(x_tilda(1,:),'LineWidth',2);
        grid on
        title('x tilda 1');
        xlabel('t');
        subplot(1,2,2)
        plot(x_tilda(2,:),'LineWidth',2);
        grid on
        title('x tilda 2');
        xlabel('t');

        figure(23);
        subplot(1,2,1)
        plot(uu_init(1,:,1),'LineWidth',2);
        grid on
        ylabel('u (Nm)');
        xlabel('t');
        title('u_1 init');
        subplot(1,2,2)
        plot(uu_init(2,:,1),'LineWidth',2);
        grid on
        ylabel('u (Nm)');
        xlabel('t');
        title('u_2 init');
    end

end