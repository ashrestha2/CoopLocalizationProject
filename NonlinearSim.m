%% ODE-45 Simulation
clc
clear
close all

% Initial conditions
L = 0.5;
deltaT = 0.1;
xi_g0 = 10;
eta_g0 = 0;
theta_g0 = pi/2;
v_g0 = 2;
phi_g0 = -pi/18;
xi_a0 = -60;
eta_a0 = 0;
theta_a0 = -pi/2;
v_a0 = 12;
w_a0 = pi/25;

x0 = [xi_g0;eta_g0;theta_g0;xi_a0;eta_a0;theta_a0];
u = [v_g0;phi_g0;v_a0;w_a0];
perturb_x0 = [0;1;0;0;0;0.1];

t_int = 0:deltaT:100;
options = odeset('RelTol',1E-12,'AbsTol',1E-12);
[t,X] = ode45(@(t,x) ugvEOM(t,x,u,L),t_int,x0,options); 

n = length(x0);
var = {'$\xi_{g}$ [m]','$\eta_{g}$ [m]','$\theta_{g}$ [rads]','$\xi_{a}$ [m]','$\eta_{a}$ [m]','$\theta_{a}$ [rads]'};
figure
for i = 1:n
    subplot(n,1,i)
    plot(t,X(:,i))
    xlabel('Time (secs)')
    ylabel(var{i},'Interpreter','latex')
end
sgtitle('States vs Time, Full Nonlinear Dynamics Simulation')