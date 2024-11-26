%% ASEN 5044 PROJECT: COOP LOCALIZATION
% Authors: Ananya Shrestha, Maggie Wussow, Cara Spencer
% Last Updated: 11/26/2024

%% Part 1: DETERMINISTIC SYSTEM ANALYSIS
%%% housekeeping
clc
clear
close all

%%% givens
L = 0.5; % length of UGV [m]
deltaT = 0.1; % time interval

% Nominal conditions
%   UGV
xi_g0 = 10; % east position [m]
eta_g0 = 0; % north position [m]
theta_g0 = pi/2; % heading angle [rad]
v_g0 = 2; % linear velocity [m/s]
phi_g0 = -pi/18; % steering angle [rad]
%   UAV
xi_a0 = -60; % east position [m]
eta_a0 = 0; % north position [m]
theta_a0 = -pi/2; % heading angle [m]
v_a0 = 12; % linear velocity [m/s]
w_a0 = pi/25; % angluar rate [rad/s]

%%% 1) CT Linearized Model

Anom = [
    0, 0, -v_g0 * sin(theta_g0), 0, 0, 0;
    0, 0,  v_g0 * cos(theta_g0), 0, 0, 0;
    0, 0,  0,                  0, 0, 0;
    0, 0,  0,                  0, -v_a0 * sin(theta_a0), 0;
    0, 0,  0,                  0,  v_a0 * cos(theta_a0), 0;
    0, 0,  0,                  0,  0, 0
];

Bnom = [
    cos(theta_g0), 0,               0, 0;
    sin(theta_g0), 0,               0, 0;
    tan(phi_g0)/L, 0, v_g0/(L*cos(phi_g0)^2), 0;
    0, cos(theta_a0),               0, 0;
    0, sin(theta_a0),               0, 0;
    0, 0,                          0, 1
];

Gammanom = [
    1, 0, 0, 0, 0, 0;
    0, 1, 0, 0, 0, 0;
    0, 0, 1, 0, 0, 0;
    0, 0, 0, 1, 0, 0;
    0, 0, 0, 0, 1, 0;
    0, 0, 0, 0, 0, 1
];

Cnom = [
    (eta_a0 - eta_g0) / ((xi_a0 - xi_g0)^2 + (eta_a0 - eta_g0)^2), ...
   -(xi_a0 - xi_g0) / ((xi_a0 - xi_g0)^2 + (eta_a0 - eta_g0)^2), 0, 0, 0, 0;
   -(eta_a0 - eta_g0) / ((xi_a0 - xi_g0)^2 + (eta_a0 - eta_g0)^2), ...
   -(xi_g0 - xi_a0) / ((xi_g0 - xi_a0)^2 + (eta_g0 - eta_a0)^2), ...
    (eta_g0 - eta_a0) / ((xi_g0 - xi_a0)^2 + (eta_g0 - eta_a0)^2), 0, 0, 0;
   -(xi_a0 - xi_g0) / ((xi_a0 - xi_g0)^2 + (eta_a0 - eta_g0)^2), ...
    (xi_g0 - xi_a0) / ((xi_g0 - xi_a0)^2 + (eta_g0 - eta_a0)^2), ...
   -(xi_g0 - xi_a0) / ((xi_g0 - xi_a0)^2 + (eta_g0 - eta_a0)^2), 0, 0, 0;
    0, 0, 0, -1, 0, 0;
    0, 0, 0,  0, 1, 0
];

Dnom = zeros(5, 6);

%%% 2) DT Linearized Model


%%% 3) Comparison of DT and Non-linear Simulations
% Initial conditions
x0 = [xi_g0;eta_g0;theta_g0;xi_a0;eta_a0;theta_a0];
u = [v_g0;phi_g0;v_a0;w_a0];
perturb_x0 = [0;1;0;0;0;0.1];

% Non-linear Simulation
t_int = 0:deltaT:100;
options = odeset('RelTol',1E-12,'AbsTol',1E-12);
[t,X] = ode45(@(t,x) ugvEOM(t,x,u,L),t_int,x0+perturb_x0,options); 

n = length(x0);
var = {'$\xi_{g}$ [m]','$\eta_{g}$ [m]','$\theta_{g}$ [rads]','$\xi_{a}$ [m]','$\eta_{a}$ [m]','$\theta_{a}$ [rads]'};
figure
for i = 1:n
    subplot(n,1,i)
    plot(t,X(:,i))
    if i == 3 || i == 6
        plot(t,wrapToPi(X(:,i)))
    end
    xlabel('Time (secs)')
    ylabel(var{i},'Interpreter','latex')
end
sgtitle('States vs Time, Full Nonlinear Dynamics Simulation')

%% PART 2: STOCHASTIC NONLINEAR FILTERING

%%% 4) Implement and tune KF

load("cooplocalization_finalproj_KFdata.mat")

% a) 