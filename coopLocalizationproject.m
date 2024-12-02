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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% 2) DT Linearized Model

% testing nominal iterations
time_steps = 0:0.1:100;
%%% trying nomial again 
theta_g = theta_g0;
theta_a = theta_a0;
nom_east_g(:,1) = xi_g0;
nom_north_g(:,1) = eta_g0;
nom_theta_g(:,1) = theta_g0;
nom_east_a(:,1) = xi_a0;
nom_north_a(:,1) = eta_a0;
nom_theta_a(:,1) = theta_a0;
for nom_step = 1:length(time_steps)-1
    nom_east_g(:,nom_step+1) = v_g0*cos(theta_g)*(deltaT*nom_step);
    nom_north_g(:,nom_step+1) = v_g0*sin(theta_g)*(deltaT*nom_step);
    nom_theta_g(:,nom_step+1) = wrapToPi( v_g0/L * tan(phi_g0) * (deltaT*nom_step)) ;
    nom_east_a(:,nom_step+1) = v_a0*cos(theta_a)*(deltaT*nom_step);
    nom_north_a(:,nom_step+1) = v_a0*sin(theta_a)*(deltaT*nom_step);
    nom_theta_a(:,nom_step+1) = wrapToPi( w_a0 * (deltaT*nom_step)) ;

    theta_a = nom_theta_a(:,nom_step);
    theta_g = nom_theta_g(:,nom_step);
end
x_nom_test = [nom_east_g, nom_north_g, nom_theta_g, nom_east_a, nom_north_a, nom_theta_a];

figure(); hold on;
plot(nom_east_g,nom_north_g,'b')
plot(nom_east_a,nom_north_a,'r')
plot(nom_east_g(:,1),nom_north_g(:,1),'go')
plot(nom_east_a(:,1),nom_north_a(:,1),'go')
xlabel('e')
ylabel('n')

figure(); hold on;
subplot(3,2,1); hold on
plot(time_steps,nom_east_g)
subplot(3,2,2); hold on
plot(time_steps,nom_east_a)
subplot(3,2,3); hold on
plot(time_steps,nom_north_g)
subplot(3,2,4); hold on
plot(time_steps,nom_north_a)
subplot(3,2,5); hold on
plot(time_steps,nom_theta_g)
subplot(3,2,6); hold on
plot(time_steps,nom_theta_a)


%%
%%% assume we have F, G, H, M 
% initial perturbations -- CHANGE LATER 
del_east_g0 = 0;
del_north_g0 = 1; 
del_theta_g0 = deg2rad(0);
del_east_a0 = 0;
del_north_a0 = 0;
del_theta_a0 = 0.1;
del_x_0 = [del_east_g0;del_north_g0; del_theta_g0; del_east_a0; del_north_a0; del_theta_a0];
%perturb_x0 = [0;1;0;0;0;0.1];

% create time steps 
end_time = 400*0.1;
time_steps = 0:0.1:end_time;

 del_x(:,1) = del_x_0;
% run through linearized perturbation iterations 
for lin_iter = 1:length(time_steps)
    % x(t) = phi(t,t0) * x(0)
    del_x(:,lin_iter+1) = F *  del_x(:,lin_iter); %assuming no control input pertubations (del_u = 0)
    % double check -- 
    del_x_check(:,lin_iter) = F^lin_iter * del_x_0; 
end

% running the nominal points
% pulled from hw 7 -- orbit given omega and del_t 
del_t = deltaT;
F_func =  @(omega) [1 (sin(omega*del_t)/omega) 0 -(1-cos(omega*del_t))/omega;
    0 cos(omega*del_t) 0 -sin(omega*del_t);
    0 (1-cos(omega*del_t))/omega 1 (sin(omega*del_t)/omega);
    0 sin(omega*del_t) 0 cos(omega*del_t)];

% find the F matricies for nominal
F_ugv = F_func(v_g0/L * tan(phi_g0)); % WILL NEED TO CHECK IF THIS IS CORRECTA for omega 
F_uav = F_func(w_a0);
% combine for a full nominal state F
F_nom = blkdiag(F_ugv,F_uav);

% find the initial nominal states
x_nom_0 = [xi_g0;eta_g0;theta_g0;xi_a0;eta_a0;theta_a0];
x_nom(:,1) = x_nom_0;

% iterate through the nominal points 
for nom_iter = 1:length(time_steps)
    x_nom(:,nom_iter+1) = F_nom * x_nom(:,nom_iter);
end
figure(); hold on;
plot(x_nom(:,1),x_nom(:,2),'b')
plot(x_nom(:,4),x_nom(:,5),'r')
xlabel('e')
ylabel('n')

%%% trying nomial again 
theta_g = theta_g0;
theta_a = theta_a0;
for nom_step = 1:length(time_steps)
    nom_east_g(:,nom_step) = v_g0*cos(theta_g)*deltaT;
    nom_north_g(:,nom_step) = v_g0*sin(theta_g)*deltaT;
    nom_theta_g(:,nom_step) = v_g0/L * tan(phi_g0) * deltaT;
    nom_east_a(:,nom_step) = v_a0*cos(theta_a)*deltaT;
    nom_north_a(:,nom_step) = v_a0*sin(theta_a)*deltaT;
    nom_theta_a(:,nom_step) = w_a0 * deltaT;

    theta_a = nom_theta_a(:,nom_step);
    theta_g = nom_theta_g(:,nom_step);
end
x_nom_test = [nom_east_g, nom_north_g, nom_theta_g, nom_east_a, nom_north_a, nom_theta_a];

% full x = nom + perturbations
x_linear = x_nom + del_x;

% plotting 
n = length(x_nom_0);
var = {'$\xi_{g}$ [m]','$\eta_{g}$ [m]','$\theta_{g}$ [rads]','$\xi_{a}$ [m]','$\eta_{a}$ [m]','$\theta_{a}$ [rads]'};
figure();
for j = 1:n
    subplot(n,1,j)
    plot(time_steps,x_linear(:,j))
    if j == 3 || j == 6
        plot(t,wrapToPi(x_linear(:,i)))
    end
    ylabel(var{j},'Interpreter','latex')
end
xlabel('Time (secs)')
sgtitle('States vs Time, Full Linearized Dynamics Simulation')
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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


% figure(); hold on;
% plot(X(:,1),X(:,2),'b')
% plot(X(:,4),X(:,5),'r')
figure(); hold on;
subplot(3,2,1); hold on
plot(t,X(:,1))
subplot(3,2,2)
plot(t,X(:,4))
subplot(3,2,3)
plot(t,X(:,2))
subplot(3,2,4)
plot(t,X(:,5))
subplot(3,2,5)
plot(t,wrapToPi(X(:,3)))
subplot(3,2,6)
plot(t,wrapToPi(X(:,6)))

%% PART 2: STOCHASTIC NONLINEAR FILTERING

%%% 4) Implement and tune KF

load("cooplocalization_finalproj_KFdata.mat")

% a) 