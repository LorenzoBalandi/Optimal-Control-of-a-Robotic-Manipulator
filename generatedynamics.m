%% Generate dynamics 
%  Using the symbolic math toolbox
%  it generates a function with its gradient and hessian
% 
% Lorenzo Sforni
% Group 18: Balandi, Ghinelli, Prandin

state_dim = 4; % 4 states x1, x2, x3, x4
input_dim = 2; % 2 input torques u1, u2

syms xx [state_dim,1] real
syms uu [input_dim,1] real
syms pp [state_dim,1] real

% parameters
mm1 = 2;
mm2 = 2;
JJ1 = 0.5;
JJ2 = 0.5;
ll1 = 1;
ll2 = 1;
rr1 = 0.5;
rr2 = 0.5;
gg = 9.81;
dt = 1e-3;
% components of the dynamics

uu = [uu(1); uu(2)]; % Input torques applied on the two joints

% Let's define the matrices M, C, g, which define the dynamic model 
MM = [mm1*rr1^2+mm2*(ll1^2+rr2^2+2*ll1*rr2*cos(xx(3)))+JJ1+JJ2, mm2*(rr2^2+ll1*rr2*cos(xx(3)))+JJ2;...
    mm2*(rr2^2+ll1*rr2*cos(xx(3)))+JJ2, mm2*rr2^2+JJ2];

CC = -mm2*ll1*rr2*sin(xx(3))*[xx(4), xx(2)+xx(4); -xx(2), 0];

GG = [(mm1*rr1+mm2*ll1)*gg*cos(xx(1))+mm2*gg*rr2*cos(xx(1)+xx(3));...
    mm2*gg*rr2*cos(xx(1)+xx(3))];


% These auxiliary variables are used in order to define the derivatives of 
% the second and fourth components of the state (theta1_dot and theta2_dot)
qq = [xx(1);xx(3)];
qq_dot = [xx(2);xx(4)];
qq_doubledot = MM\(uu - CC*qq_dot - GG);

xx_1_dot = xx(2);
xx_2_dot = qq_doubledot(1);
xx_3_dot = xx(4);
xx_4_dot = qq_doubledot(2);

% Discretization through Euler method
f1 = xx(1) + dt * xx_1_dot;
f2 = xx(2) + dt * xx_2_dot;
f3 = xx(3) + dt * xx_3_dot;
f4 = xx(4) + dt * xx_4_dot;

% xxp xontains the values of the state x at time t+1
xxp = [f1; f2; f3; f4];

dfx = jacobian(xxp,xx)';
dfu = jacobian(xxp,uu)';

% Generate the hessian matrices
d2fxx = cell(1, state_dim);
d2fuu = cell(1, state_dim);
d2fux = cell(1, state_dim);

pfxx = zeros(state_dim,state_dim);
pfuu = zeros(input_dim,input_dim);
pfux = zeros(input_dim,state_dim)';

for ii = 1:state_dim
    d2fxx{ii} = jacobian(dfx(:,ii), xx);
    pfxx = pfxx + (d2fxx{ii}*pp(ii));
    
    d2fuu{ii} = jacobian(dfu(:,ii), uu);
    pfuu = pfuu + (d2fuu{ii}*pp(ii));
    
    d2fux{ii} = jacobian(dfx(:,ii), uu);
    pfux = pfux + (d2fux{ii}*pp(ii));
end
pfux=pfux';

matlabFunction(xxp, dfx, dfu, pfxx, pfuu, pfux, 'File', 'dynamics','Vars',{xx,uu,pp});