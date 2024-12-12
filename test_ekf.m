clear
clc
close all

load("cooplocalization_finalproj_KFdata.mat")

% State transition function (f)
f = @(x, u, dt) [
    x(1) + u(1)*cos(x(3))*dt;
    x(2) + u(1)*sin(x(3))*dt;
    x(3) + (u(1)/0.5)*tan(u(2))*dt;
    x(4) + u(3)*cos(x(6))*dt;
    x(5) + u(3)*sin(x(6))*dt;
    x(6) + u(4)*dt;
];

% Measurement function (h)
h = @(x) [
    atan2(x(5) - x(2), x(4) - x(1)) - x(3);
    sqrt((x(4) - x(1))^2 + (x(5) - x(2))^2);
    atan2(x(2) - x(5), x(1) - x(4)) - x(6);
    x(4);
    x(5);
];

% State transition Jacobian (F_func)
F_func = @(x, u, dt) [
    1, 0, -u(1)*sin(x(3))*dt, 0, 0, 0;
    0, 1,  u(1)*cos(x(3))*dt, 0, 0, 0;
    0, 0,  1,                 0, 0, 0;
    0, 0,  0,                 1, 0, -u(3)*sin(x(6))*dt;
    0, 0,  0,                 0, 1,  u(3)*cos(x(6))*dt;
    0, 0,  0,                 0, 0,  1
];

% Measurement Jacobian (H_func)
H_func = @(x) [
    (x(5) - x(2))/((x(4) - x(1))^2 + (x(5) - x(2))^2),     -(x(4) - x(1))/((x(4) - x(1))^2 + (x(5) - x(2))^2),     -1, -(x(5) - x(2))/((x(4) - x(1))^2 + (x(5) - x(2))^2),     (x(4) - x(1))/((x(4) - x(1))^2 + (x(5) - x(2))^2),      0;
    (x(1) - x(4))/sqrt((x(1) - x(4))^2 + (x(2) - x(5))^2),  (x(2) - x(5))/sqrt((x(1) - x(4))^2 + (x(2) - x(5))^2),   0, -(x(1) - x(4))/sqrt((x(1) - x(4))^2 + (x(2) - x(5))^2), -(x(2) - x(5))/sqrt((x(1) - x(4))^2 + (x(2) - x(5))^2), 0;
    -(x(2) - x(5))/((x(1) - x(4))^2 + (x(2) - x(5))^2),     (x(1) - x(4))/((x(1) - x(4))^2 + (x(2) - x(5))^2),       0, (x(2) - x(5))/((x(1) - x(4))^2 + (x(2) - x(5))^2),      -(x(1) - x(4))/((x(1) - x(4))^2 + (x(2) - x(5))^2),     -1;
    0, 0, 0, 1, 0, 0;
    0, 0, 0, 0, 1, 0
    ];

% make a struct for constants
const.L = 0.5; % length of UGV [m]
const.deltaT = 0.1; % time interval
%UGV
const.v_g0 = 2; % linear velocity [m/s]
const.phi_g0 = -pi/18; % steering angle [rad]
%UAV
const.v_a0 = 12; % linear velocity [m/s]
const.w_a0 = pi/25; % angluar rate [rad/s]

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
x0 = [xi_g0;eta_g0;theta_g0;xi_a0;eta_a0;theta_a0];
const.x0 = x0;
endTime = 100;
P0 = diag([5, 5, pi/4, 10, 10, pi/4]);

rng(100)

[time_iter,x_noisy, y_noisy] = TMTSim(const, Qtrue,Rtrue,endTime);

% Define parameters
y_meas = y_noisy; % Measurements (5 x 100)
x0 = [10; 0; pi/2; -60; 0; -pi/2]; % Initial state
Q = 10 * eye(6); % Process noise covariance
R = Rtrue; % Measurement noise covariance
u = repmat([2; -pi/18; 12; pi/25], 1, 1000); % Constant control inputs
dt = 0.1; % Sampling time

N = 1000; % Number of time steps
t = 0:dt:endTime;

% Call EKF
[x_plus, P_plus, sigma, innovation, Sk, y_calc, F_matrices] = ekf(y_meas, x0, P0, Q, R, u, dt, N, f, h, F_func, H_func);

n = length(x0);
var = {'$\xi_{g}$ [m]','$\eta_{g}$ [m]','$\theta_{g}$ [rads]','$\xi_{a}$ [m]','$\eta_{a}$ [m]','$\theta_{a}$ [rads]'};
figure
for i = 1:n
    subplot(n,1,i)
    hold on
    if i == 3 || i == 6
        plot(t,wrapToPi(x_plus(i,:)),'r',LineWidth=1.2)
        plot(t,wrapToPi(x_noisy(:,i)),'b',LineWidth=1.2)
    else
        plot(t,x_plus(i,:),'r',LineWidth=1.2)
        plot(t,x_noisy(:,i),'b',LineWidth=1.2)
    end
    ylabel(var{i},'Interpreter','latex')
end
xlabel('Time (secs)','Interpreter','latex')
sgtitle('States vs Time, Full Nonlinear Dynamics Simulation','Interpreter','latex')

p = min(size(y_noisy));
var = {'$\gamma_{ag}$ [rads]','$\rho_{g}$ [m]','$\gamma_{ga}$ [rads]','$\xi_{a}$ [m]','$\eta_{a}$ [m]'};
figure();

for i = 1:p
    subplot(p,1,i)
    hold on
    if i == 1 || i == 3
        plot(t(2:end),wrapToPi(y_calc(i,:)),'r',LineWidth=1.2)
        plot(t(2:end),wrapToPi(y_noisy(i,:)),'b',LineWidth=1.2)
    else
        plot(t(2:end),y_calc(i,:),'r',LineWidth=1.2)
        plot(t(2:end),y_noisy(i,:),'b',LineWidth=1.2)
    end
    ylabel(var{i},'Interpreter','latex')
end
xlabel('Time (secs)','Interpreter','latex')
sgtitle('Measurements vs Time, Full Nonlinear Dynamics Simulation','Interpreter','latex')

%plot errors
error_x = x_plus - x_noisy';
error_y = y_calc - y_noisy;

% plotting error with cov bounds 
var = {'$e_{\xi_{g}}$ [m]','$e_{\eta_{g}}$ [m]','$e_{\theta_{g}}$ [rads]','$e_{\xi_{a}}$ [m]','$e_{\eta_{a}}$ [m]','$e_{\theta_{a}}$ [rads]'};
figure(20);
for i = 1:n
    subplot(n,1,i); hold on;
    if i == 3 || i == 6
        plot(wrapToPi(error_x(i,:)),'r',LineWidth=1.2)
        plot(wrapToPi(2*sigma(i,:)),'b--',LineWidth=1.2)
        plot(wrapToPi(-2*sigma(i,:)),'b--',LineWidth=1.2)
    else
        plot(error_x(i,:),'r',LineWidth=1.2)
        plot(2*sigma(i,:),'b--',LineWidth=1.2)
        plot(-2*sigma(i,:),'b--',LineWidth=1.2)
    end
    ylabel(var{i},'Interpreter','latex')
end
xlabel('Time (secs)','Interpreter','latex')
sgtitle('States Error vs Time','Interpreter','latex')

var = {'$e_{\gamma_{ag}}$ [rads]','$e_{\rho_{g}}$ [m]','$e_{\gamma_{ga}}$ [rads]','$e_{\xi_{a}}$ [m]','$e_{\eta_{a}}$ [m]'};
figure(21);
for i = 1:p
    subplot(p,1,i); hold on;
    if i == 1 || i == 3
        plot(wrapToPi(innovation(i,:)),'g',LineWidth=1.2)
        plot(wrapToPi(-2*sqrt(reshape(Sk(i,i,2:end),[1,length(Sk)-1]))),'b--',LineWidth=1.2)
        plot(wrapToPi(+2*sqrt(reshape(Sk(i,i,2:end),[1,length(Sk)-1]))),'b--',LineWidth=1.2)
    else
        plot(innovation(i,:),'g',LineWidth=1.2)
        plot(-2*sqrt(reshape(Sk(i,i,2:end),[1,length(Sk)-1])),'b--',LineWidth=1.2)
        plot(+2*sqrt(reshape(Sk(i,i,2:end),[1,length(Sk)-1])),'b--',LineWidth=1.2)  
    end
    ylabel(var{i},'Interpreter','latex')
end
xlabel('Time (secs)','Interpreter','latex')
sgtitle('Measurements Error vs Time','Interpreter','latex')
%% MOnte Carlo RUn
num = 5;

for i = 1:num
    % Generate true trajectory and measurements from system
    [t,xtrue,ytrue] = TMTSim(const,Qtrue,Rtrue,endTime);

    [x_plus, P_plus, sigma, innovation, Sk, y_calc, F_matrices] = ekf(y_meas, x0, P0, Q, R, u, dt, N, f, h, F_func, H_func);
    
    % xtrue(:,3) = wrapToPi(xtrue(:,3));
    % xtrue(:,6) = wrapToPi(xtrue(:,6));
    % x_plus(3,:) = wrapToPi(x_plus(3,:));
    % x_plus(6,:) = wrapToPi(x_plus(6,:));

    % Calculate NEES and NIS
    for j = 1: length(t)
        error_x(:,i,j) = xtrue(j,:) - x_plus(:,j)';
        error_x(3,i,j) = wrapToPi(error_x(3,i,j));
        error_x(6,i,j) = wrapToPi(error_x(6,i,j));
        NEES(i,j) = (error_x(:,i,j)')*inv(P_plus(:,:,j))*(error_x(:,i,j));
        if j > 2
            NIS(i,j-1) = (innovation(:,j-1))'*inv(Sk(:,:,j-1))*(innovation(:,j-1));
        end
    end
end

% NEES Test:
epsNEESbar = mean(NEES,1);
alphaNEES = 0.01;
Nnx = num*n;
r1x = chi2inv(alphaNEES/2, Nnx)./num;
r2x = chi2inv(1-alphaNEES/2, Nnx)./num;

figure
plot(epsNEESbar,'ro','MarkerSize',6,'LineWidth',2),hold on
plot(r1x*ones(size(epsNEESbar)),'r--','LineWidth',2)
plot(r2x*ones(size(epsNEESbar)),'r--','LineWidth',2)
ylabel('NEES statistic, $\bar{\epsilon}_x$','FontSize',14,'Interpreter','latex')
xlabel('time step, k','FontSize',14)
title('NEES Estimation Results','FontSize',14)
legend('NEES @ time k', 'r_1 bound', 'r_2 bound')
ylim([0 15])

epsNISbar = mean(NIS,1);
alphaNIS = 0.1;
Nny = num*p;
r1y = chi2inv(alphaNIS/2,Nny)./num;
r2y = chi2inv(1-alphaNIS/2,Nny)./num;

figure
plot(epsNISbar,'bo','MarkerSize',6,'LineWidth',2),hold on
plot(r1y*ones(size(epsNISbar)),'b--','LineWidth',2)
plot(r2y*ones(size(epsNISbar)),'b--','LineWidth',2)
ylabel('NIS statistic, $\bar{\epsilon}_y$','FontSize',14,'Interpreter','latex')
xlabel('time step, k','FontSize',14)
title('NIS Estimation Results','FontSize',14)
legend('NIS @ time k', 'r_1 bound', 'r_2 bound')
ylim([0 15])

% % plot innovation
% figure(17);
% subplot(5,1,1); hold on;
% plot(innovation(1,:));
% plot(+2*sqrt(reshape(Sk(1,1,2:end),[1,length(Sk)-1])),'b--',LineWidth=1.2)
% plot(-2*sqrt(reshape(Sk(1,1,2:end),[1,length(Sk)-1])),'b--',LineWidth=1.2)
% subplot(5,1,2); hold on;
% plot(innovation(2,:));
% plot(+2*sqrt(reshape(Sk(2,2,2:end),[1,length(Sk)-1])),'b--',LineWidth=1.2)
% plot(-2*sqrt(reshape(Sk(2,2,2:end),[1,length(Sk)-1])),'b--',LineWidth=1.2)
% subplot(5,1,3); hold on;
% plot(innovation(3,:));
% plot(+2*sqrt(reshape(Sk(3,3,2:end),[1,length(Sk)-1])),'b--',LineWidth=1.2)
% plot(-2*sqrt(reshape(Sk(3,3,2:end),[1,length(Sk)-1])),'b--',LineWidth=1.2)
% subplot(5,1,4); hold on;
% plot(innovation(4,:));
% plot(+2*sqrt(reshape(Sk(4,4,2:end),[1,length(Sk)-1])),'b--',LineWidth=1.2)
% plot(-2*sqrt(reshape(Sk(4,4,2:end),[1,length(Sk)-1])),'b--',LineWidth=1.2)
% subplot(5,1,5); hold on;
% plot(innovation(5,:));
% sgtitle('innovation')


