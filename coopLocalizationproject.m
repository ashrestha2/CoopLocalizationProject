%% ASEN 5044 PROJECT: COOP LOCALIZATION
% Authors: Ananya Shrestha, Maggie Wussow, Cara Spencer
% Last Updated: 12/7/2024

%% Part 1: DETERMINISTIC SYSTEM ANALYSIS
%%% housekeeping
clc
clear
close all

%%% givens
L = 0.5; % length of UGV [m]
deltaT = 0.1; % time interval

% Initial nominal conditions
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1) CT Linearized Model
syms eta_g xi_g theta_g v_g phi_g eta_a xi_a theta_a v_a w_a
Anom = [
    0, 0, -v_g * sin(theta_g), 0, 0, 0;
    0, 0,  v_g * cos(theta_g), 0, 0, 0;
    0, 0,  0,                  0, 0, 0;
    0, 0,  0,                  0, 0, -v_a * sin(theta_a);
    0, 0,  0,                  0, 0, v_a * cos(theta_a);
    0, 0,  0,                  0,  0, 0
];

Bnom = [
    cos(theta_g), 0, 0, 0;
    sin(theta_g), 0, 0, 0;
    tan(phi_g)/L, 0, v_g/(L*sec(phi_g)^2), 0;
    0, 0, cos(theta_a), 0;
    0, 0, sin(theta_a), 0;
    0, 0, 0, 1
];

C = [ 
    (eta_a - eta_g) / ((xi_a - xi_g)^2 + (eta_a - eta_g)^2), ...
    -(xi_a - xi_g) / ((xi_a - xi_g)^2 + (eta_a - eta_g)^2), ...
    -1, ...
    -(eta_a - eta_g) / ((xi_a - xi_g)^2 + (eta_a - eta_g)^2), ...
    (xi_a - xi_g) / ((xi_a - xi_g)^2 + (eta_a - eta_g)^2), ...
    0;

    (xi_g - xi_a) / ((xi_g - xi_a)^2 + (eta_g - eta_a)^2), ...
    (eta_g - eta_a) / ((xi_g - xi_a)^2 + (eta_g - eta_a)^2), ...
    0, ...
    -(xi_g - xi_a) / ((xi_g - xi_a)^2 + (eta_g - eta_a)^2), ...
    (xi_g - xi_a) / ((xi_g - xi_a)^2 + (eta_g - eta_a)^2), ...
    0;

    -(eta_g - eta_a) / ((xi_g - xi_a)^2 + (eta_g - eta_a)^2), ...
    (xi_g - xi_a) / ((xi_g - xi_a)^2 + (eta_g - eta_a)^2), ...
    -(xi_g - xi_a) / ((xi_g - xi_a)^2 + (eta_g - eta_a)^2), ...
    (eta_g - eta_a) / ((xi_g - xi_a)^2 + (eta_g - eta_a)^2), ...
    (xi_g - xi_a) / ((xi_g - xi_a)^2 + (eta_g - eta_a)^2), ...
    -1;

    0, 0, 0, 1, 0, 0;

    0, 0, 0, 0, 1, 0
];

Gammanom = [
    1, 0, 0, 0, 0, 0;
    0, 1, 0, 0, 0, 0;
    0, 0, 1, 0, 0, 0;
    0, 0, 0, 1, 0, 0;
    0, 0, 0, 0, 1, 0;
    0, 0, 0, 0, 0, 1
];

Dnom = zeros(5,6);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2) DT Linearized Model

%%%% Finding the Nominal trajectory 
% set the time steps
time_steps = 0:deltaT:100;
% initial state
x0 = [xi_g0;eta_g0;theta_g0;xi_a0;eta_a0;theta_a0];
% run ode45 to find the nominal trajectory at every time step
[t_nom,x_nom] = ode45(@(t,y) FindNominal(t, y, v_g0, v_a0, L, phi_g0, w_a0),time_steps,x0,[]);

% % plotting jsut to check
% figure(); hold on;
% subplot(3,2,1); hold on
% plot(time_steps,x_nom(:,1))
% subplot(3,2,2); hold on
% plot(time_steps,x_nom(:,4))
% subplot(3,2,3); hold on
% plot(time_steps,x_nom(:,2))
% subplot(3,2,4); hold on
% plot(time_steps,x_nom(:,5))
% subplot(3,2,5); hold on
% plot(time_steps,wrapToPi(x_nom(:,3)))
% subplot(3,2,6); hold on
% plot(time_steps,wrapToPi(x_nom(:,6)))
% title('Nominal Trajectory')

%initial perturbation
perturb_x0 = [0;1;0;0;0;0.1];
del_x(:,1) = perturb_x0;

% run through linearized perturbation iterations 
for lin_iter = 1:length(time_steps)-1
    % pull the F,G,H,M matricies solved at a different nominal point
    [F,G,H,M] = CT_to_DT(x_nom(lin_iter,:),L,v_g0,v_a0,phi_g0,w_a0,deltaT);

    % x(t) = phi(t,t0) * x(0)
    del_x(:,lin_iter+1) = F *  del_x(:,lin_iter); %assuming no control input pertubations (del_u = 0)
    % double check -- 
    del_x_check(:,lin_iter) = F^lin_iter * perturb_x0; 

    % find the linearlized measurements
    del_y(:,lin_iter) = H * del_x(:,lin_iter);
end

% full x = nom + perturbations
x_linear = x_nom + del_x';

% plotting 
n = length(x0);
var = {'$\xi_{g}$ [m]','$\eta_{g}$ [m]','$\theta_{g}$ [rads]','$\xi_{a}$ [m]','$\eta_{a}$ [m]','$\theta_{a}$ [rads]'};
figure();
for j = 1:n
    subplot(n,1,j)
    plot(time_steps,x_linear(:,j),'b',LineWidth=1.2)
    if j == 3 || j == 6
        plot(time_steps,wrapToPi(x_linear(:,j)),'b',LineWidth=1.2)
    end
    ylabel(var{j},'Interpreter','latex')
end
xlabel('Time (secs)','Interpreter','latex')
sgtitle('States vs Time, Full Linearized Dynamics Simulation','Interpreter','latex')

% find the nominal measurements 
y_nom = findYnom(x_nom(2:end,:));

% find the total measurements 
y_linear = y_nom + del_y;

p = min(size(y_nom));
var = {'$\gamma_{ag}$ [rads]','$\rho_{g}$ [m]','$\gamma_{ga}$ [rads]','$\xi_{a}$ [m]','$\eta_{a}$ [m]'};
figure();
for j = 1:p
    subplot(p,1,j)
    plot(time_steps(2:end),y_linear(j,:),'b',LineWidth=1.2)
    if j == 1 || j == 3
        plot(time_steps(2:end),wrapToPi(y_linear(j,:)),'b',LineWidth=1.2)
    end
    ylabel(var{j},'Interpreter','latex')
end
xlabel('Time (secs)','Interpreter','latex')
sgtitle('Measurements vs Time, Full Linearized Dynamics Simulation','Interpreter','latex')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3) Comparison of DT and Non-linear Simulations
% Initial conditions
x0 = [xi_g0;eta_g0;theta_g0;xi_a0;eta_a0;theta_a0];
u = [v_g0;phi_g0;v_a0;w_a0];
w = [0;0;0;0;0;0];
perturb_x0 = [0;1;0;0;0;0.1];

% Non-linear Simulation
t_int = 0:deltaT:100;
options = odeset('RelTol',1E-12,'AbsTol',1E-12);
[t,X] = ode45(@(t,x) ugvEOM(t,x,u,w,L),t_int,x0+perturb_x0,options); 

% Non-linear Measurement Simulation
Y = findYnom(X(2:end,:));

n = length(x0);
var = {'$\xi_{g}$ [m]','$\eta_{g}$ [m]','$\theta_{g}$ [rads]','$\xi_{a}$ [m]','$\eta_{a}$ [m]','$\theta_{a}$ [rads]'};
figure
for i = 1:n
    subplot(n,1,i)
    plot(t,X(:,i),'r',LineWidth=1.2)
    if i == 3 || i == 6
        plot(t,wrapToPi(X(:,i)),'r',LineWidth=1.2)
    end
    ylabel(var{i},'Interpreter','latex')
end
xlabel('Time (secs)','Interpreter','latex')
sgtitle('States vs Time, Full Nonlinear Dynamics Simulation','Interpreter','latex')

p = min(size(y_nom));
var = {'$\gamma_{ag}$ [rads]','$\rho_{g}$ [m]','$\gamma_{ga}$ [rads]','$\xi_{a}$ [m]','$\eta_{a}$ [m]'};
figure();
for i = 1:p
    subplot(p,1,i)
    plot(t(2:end),Y(i,:),'r',LineWidth=1.2)
    if i == 1 || i == 3
        plot(t(2:end),wrapToPi(Y(i,:)),'r',LineWidth=1.2)
    end
    ylabel(var{i},'Interpreter','latex')
end
xlabel('Time (secs)','Interpreter','latex')
sgtitle('Measurements vs Time, Full Nonlinear Dynamics Simulation','Interpreter','latex')

% Calculating error
error_x = x_linear - X;
error_y = y_linear - Y;

n = length(x0);
var = {'$e_{\xi_{g}}$ [m]','$e_{\eta_{g}}$ [m]','$e_{\theta_{g}}$ [rads]','$e_{\xi_{a}}$ [m]','$e_{\eta_{a}}$ [m]','$e_{\theta_{a}}$ [rads]'};
figure
for i = 1:n
    subplot(n,1,i)
    plot(t,error_x(:,i),'g',LineWidth=1.2)
    if i == 3 || i == 6
        plot(t,wrapToPi(error_x(:,i)),'g',LineWidth=1.2)
    end
    ylabel(var{i},'Interpreter','latex')
end
xlabel('Time (secs)','Interpreter','latex')
sgtitle('States Error vs Time','Interpreter','latex')

p = min(size(y_nom));
var = {'$e_{\gamma_{ag}}$ [rads]','$e_{\rho_{g}}$ [m]','$e_{\gamma_{ga}}$ [rads]','$e_{\xi_{a}}$ [m]','$e_{\eta_{a}}$ [m]'};
figure();
for i = 1:p
    subplot(p,1,i)
    plot(t(2:end),error_y(i,:),'g',LineWidth=1.2)
    if i == 1 || i == 3
        plot(t(2:end),wrapToPi(error_y(i,:)),'g',LineWidth=1.2)
    end
    ylabel(var{i},'Interpreter','latex')
end
xlabel('Time (secs)','Interpreter','latex')
sgtitle('Measurements Error vs Time','Interpreter','latex')

%% PART 2: STOCHASTIC NONLINEAR FILTERING

close all;

%%% 4) Implement and tune KF

load("cooplocalization_finalproj_KFdata.mat")

% make a struct for constants
const.L = 0.5; % length of UGV [m]
const.deltaT = 0.1; % time interval
%UGV
const.v_g0 = 2; % linear velocity [m/s]
const.phi_g0 = -pi/18; % steering angle [rad]
%UAV
const.v_a0 = 12; % linear velocity [m/s]
const.w_a0 = pi/25; % angluar rate [rad/s]

const.x0 = x0;
% a) 

%%% LKF 
% IC 
del_x0 = [0;1;0;0;0;0.1];
P0 = 100 * eye(length(del_x0));
%P0 = diag([100, 100, 2*pi, 100, 100, 2*pi]);
% y_nom = findYnom(x_nom);
T = length(ydata);
LKF_time = 0:1000;
N = 5;

[epsNEESbar,r1x,r2x,epsNISbar,r1y,r2y, NEES, NIS] = FindNISNESS(N,del_x0,P0,x_nom,y_nom,@CT_to_DT,const,Qtrue,Rtrue);

% get noisy measurements and proces s
[time_iter,x_noisy, y_noisy] = TMTSim(const, Qtrue,Rtrue);

% % run batch ls -- to get initial values for LKF
% [x_bls, P_bls] = BatchLS(y_noisy,Rtrue,x_noisy,const,@CT_to_DT);

% run LKF
%[x_LKF_full, P_plus, innovation, y_LKF_total] = LKF(del_x0, P0, const, @CT_to_DT, x_nom, y_nom, y_noisy, Qtrue, Rtrue);
[x_LKF_full, P_plus, innovation, y_LKF_total] = LKF(del_x0, P0, const, @CT_to_DT, x_nom, y_nom, y_noisy, Qtrue, Rtrue);


% plotting 
n = length(del_x0);
var = {'$\xi_{g}$ [m]','$\eta_{g}$ [m]','$\theta_{g}$ [rads]','$\xi_{a}$ [m]','$\eta_{a}$ [m]','$\theta_{a}$ [rads]'};
figure();
for j = 1:n
    subplot(n,1,j); hold on;
    if j == 3 || j == 6
        plot(LKF_time,wrapToPi(x_LKF_full(j,:)),'b',LineWidth=1.2)
        plot(LKF_time,wrapToPi(x_noisy(:,j)),'r--',LineWidth=1.2)
    else
        plot(LKF_time,x_LKF_full(j,:),'b',LineWidth=1.2)
        plot(LKF_time,x_noisy(:,j),'r--',LineWidth=1.2)
    end
    ylabel(var{j},'Interpreter','latex')
end
xlabel('Time (secs)','Interpreter','latex')
sgtitle('States vs Time, LKF Simulation','Interpreter','latex')

p = min(size(y_nom));
var = {'$\gamma_{ag}$ [rads]','$\rho_{g}$ [m]','$\gamma_{ga}$ [rads]','$\xi_{a}$ [m]','$\eta_{a}$ [m]'};
figure();
for j = 1:p
    subplot(p,1,j); hold on;
    if j == 1 || j == 3
        plot(LKF_time(2:end),wrapToPi(y_LKF_total(j,:)),'b',LineWidth=1.2)
        plot(LKF_time(2:end),wrapToPi(y_noisy(j,:)),'r--',LineWidth=1.2)
    else
        plot(LKF_time(2:end),y_LKF_total(j,:),'b',LineWidth=1.2)
        plot(LKF_time(2:end),y_noisy(j,:),'r--',LineWidth=1.2)
    end
    ylabel(var{j},'Interpreter','latex')
end
xlabel('Time (secs)','Interpreter','latex')
sgtitle('Measurements vs Time, LKF Simulation','Interpreter','latex')
const.x0 = x0;

% Calculating error
error_x = x_LKF_full' - x_noisy;
error_y = y_LKF_total - y_noisy;

n = length(x0);
var = {'$e_{\xi_{g}}$ [m]','$e_{\eta_{g}}$ [m]','$e_{\theta_{g}}$ [rads]','$e_{\xi_{a}}$ [m]','$e_{\eta_{a}}$ [m]','$e_{\theta_{a}}$ [rads]'};
figure
for i = 1:n
    subplot(n,1,i)
    plot(t,error_x(:,i),'g',LineWidth=1.2)
    if i == 3 || i == 6
        plot(t,wrapToPi(error_x(:,i)),'g',LineWidth=1.2)
    end
    ylabel(var{i},'Interpreter','latex')
end
xlabel('Time (secs)','Interpreter','latex')
sgtitle('States Error vs Time','Interpreter','latex')

p = min(size(y_nom));
var = {'$e_{\gamma_{ag}}$ [rads]','$e_{\rho_{g}}$ [m]','$e_{\gamma_{ga}}$ [rads]','$e_{\xi_{a}}$ [m]','$e_{\eta_{a}}$ [m]'};
figure();
for i = 1:p
    subplot(p,1,i)
    plot(t(2:end),error_y(i,:),'g',LineWidth=1.2)
    if i == 1 || i == 3
        plot(t(2:end),wrapToPi(error_y(i,:)),'g',LineWidth=1.2)
    end
    ylabel(var{i},'Interpreter','latex')
end
xlabel('Time (secs)','Interpreter','latex')
sgtitle('Measurements Error vs Time','Interpreter','latex')


% testing NEES & NIS
for ts = 1:length(t)
    NEES = error_x(ts,:)  * inv(P_plus(:,:,ts)) * error_x(ts,:)';
    if (ts > 2)
        NIS = error_y(:,ts-1)  * inv(innovation(:,:,ts-1)) * error_y(:,ts-1)';
    end 
end


% b) 

%%% Verifying TMT
% [tmont,Xmont,Ymont] = TMTSim(const,Qtrue,Rtrue);
% 
% var = {'$\xi_{g}$ [m]','$\eta_{g}$ [m]','$\theta_{g}$ [rads]','$\xi_{a}$ [m]','$\eta_{a}$ [m]','$\theta_{a}$ [rads]'};
% figure
% for i = 1:n
%     subplot(n,1,i)
%     plot(tmont,Xmont(:,i),'r',LineWidth=1.2)
%     if i == 3 || i == 6
%         plot(tmont,wrapToPi(Xmont(:,i)),'r',LineWidth=1.2)
%     end
%     ylabel(var{i},'Interpreter','latex')
% end
% xlabel('Time (secs)','Interpreter','latex')
% sgtitle('States vs Time, Monte Carlo Simulation','Interpreter','latex')
% 
% p = min(size(y_nom));
% var = {'$\gamma_{ag}$ [rads]','$\rho_{g}$ [m]','$\gamma_{ga}$ [rads]','$\xi_{a}$ [m]','$\eta_{a}$ [m]'};
% figure();
% for i = 1:p
%     subplot(p,1,i)
%     plot(tmont(2:end),Ymont(i,:),'r',LineWidth=1.2)
%     if i == 1 || i == 3
%         plot(tmont(2:end),wrapToPi(Ymont(i,:)),'r',LineWidth=1.2)
%     end
%     ylabel(var{i},'Interpreter','latex')
% end
% xlabel('Time (secs)','Interpreter','latex')
% sgtitle('Measurements vs Time, Monte Carlo Simulation','Interpreter','latex')



