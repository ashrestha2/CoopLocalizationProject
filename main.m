%% COOP-LOCALIZATION PROJECT
clc
clear
close all

%%% GIVENS:
L = 0.5; 
v_g = 2; 
phi_g = -pi/18;
v_a = 12;
omega_a = pi/25;
deltaT = 0.1;

load('cooplocalization_finalproj_KFdata.mat')

const.L = L;
const.v_g = v_g;
const.phi_g = phi_g;
const.v_a = v_a;
const.omega_a = omega_a;
const.deltaT = deltaT;
const.endTime = 100;

xi_g0 = 10;
eta_g0 = 0;
theta_g0 = pi/2;
xi_a0 = -60;
eta_a0 = 0;
theta_a0 = -pi/2;

t = 0:deltaT:const.endTime;

xnom0 = [xi_g0;eta_g0;theta_g0;xi_a0;eta_a0;theta_a0];
delta_x0 = [0;1;0;0;0;0.1];
x0 = xnom0 + delta_x0;
xnoise(:,1) = x0;
n = length(xnom0);

h = @(x) [
    atan2(x(5) - x(2), x(4) - x(1)) - x(3);
    sqrt((x(4) - x(1))^2 + (x(5) - x(2))^2);
    atan2(x(2) - x(5), x(1) - x(4)) - x(6);
    x(4);
    x(5);
];

rng(100)

% Nominal trajectories (w/o noise or perturbations)
w = zeros(6,1);
options = odeset('RelTol',1E-12,'AbsTol',1E-12);
[tnom,xnom] = ode45(@(t,x) dubinsEOM(t,x,w,const),t,xnom0,options);
ynom = calcY(xnom(2:end,:));

% Nonlinear Dynamics with perturbation
options = odeset('RelTol',1E-12,'AbsTol',1E-12);
[tnonlin,xnonlin] = ode45(@(t,x) dubinsEOM(t,x,w,const),t,x0,options);
ynonlin = calcY(xnonlin(2:end,:));
p = min(size(ynom,1));

for i = 1: length(t)-1
    % Nonlinear Dynamics with perturbation + noise
    dt = [t(i) t(i+1)];
    options = odeset('RelTol',1E-12,'AbsTol',1E-12);
    w = mvnrnd(zeros(n,1),Qtrue)';
    [tnoise,Xnoise] = ode45(@(t,x) dubinsEOM(t,x,w,const),dt,xnoise(:,i),options);
    xnoise(:,i+1) = Xnoise(end,:)';
    v = mvnrnd(zeros(p,1),Rtrue)';
    ynoise(:,i) = h(xnoise(:,i+1)) + v;
end

% Tuning 
P0 = diag([5, 5, pi/12, 10, 10, pi/12]);
Q = 80*eye(6);

[x_plus, P_plus, sigma, y_plus, delta_y_minus, S, innov_cov] = LFK(delta_x0,P0,const,Q,Rtrue,xnom,ynom,ynoise);

P0e = 4*diag([5, 5, pi/12, 10, 10, pi/12]);
% Initialize Qekf to zeros and explicitly adjust the first and second diagonal elements
Qekf = zeros(size(P0e));
Qekf(1,1) = 10 * P0e(1,1); % Adjust first condition
Qekf(2,2) = 0.05 * P0e(2,2); % Adjust second condition
Qekf(4,4) = 0.01 * P0e(4,4); % Adjust first condition
Qekf(5,5) = 0.01 * P0e(5,5); % Adjust second condition

[x_pluse, P_plus, sigmae, y_pluse, innovation, S, innov_cove] = EFK(x0,P0e,const,Qekf,Rtrue,ynoise,h);

error_x = x_plus - xnoise;
error_xe = x_pluse - xnoise;

var = {'$\xi_{g}$ [m]','$\eta_{g}$ [m]','$\theta_{g}$ [rads]','$\xi_{a}$ [m]','$\eta_{a}$ [m]','$\theta_{a}$ [rads]'};
figure
for i = 1:n
    subplot(n,1,i)
    hold on
    if i == 3 || i == 6
        plot(t,wrapToPi(xnom(:,i)),'r',LineWidth=1.2)
        plot(t,wrapToPi(xnonlin(:,i)),'b',LineWidth=1.2)
        plot(t,wrapToPi(xnoise(i,:)),'m',LineWidth=1.2)
        plot(t,wrapToPi(x_plus(i,:)),'g',LineWidth=1.2)
        plot(t,wrapToPi(x_pluse(i,:)),'k',LineWidth=1.2)
    else
        plot(t,xnom(:,i),'r',LineWidth=1.2)
        plot(t,xnonlin(:,i),'b',LineWidth=1.2)
        plot(t,xnoise(i,:),'m',LineWidth=1.2)
        plot(t,x_plus(i,:),'g',LineWidth=1.2)
        plot(t,x_pluse(i,:),'k',LineWidth=1.2)
    end
    ylabel(var{i},'Interpreter','latex')
    legend('Nominal', 'Perturbed','Noise','LKF','EKF')
end
xlabel('Time (secs)','Interpreter','latex')
sgtitle('States vs Time','Interpreter','latex')

var = {'$\gamma_{ag}$ [rads]','$\rho_{g}$ [m]','$\gamma_{ga}$ [rads]','$\xi_{a}$ [m]','$\eta_{a}$ [m]'};
figure
for i = 1:p
    subplot(p,1,i)
    hold on
    if i == 1 || i == 3
        plot(t(2:end),wrapToPi(ynom(i,:)),'r',LineWidth=1.2)
        plot(t(2:end),wrapToPi(ynonlin(i,:)),'b',LineWidth=1.2)
        plot(t(2:end),wrapToPi(ynoise(i,:)),'m',LineWidth=1.2)
        plot(t(2:end),wrapToPi(y_plus(i,:)),'g',LineWidth=1.2)
        plot(t(2:end),wrapToPi(y_pluse(i,:)),'k',LineWidth=1.2)
    else
        plot(t(2:end),ynom(i,:),'r',LineWidth=1.2)
        plot(t(2:end),ynonlin(i,:),'b',LineWidth=1.2)
        plot(t(2:end),ynoise(i,:),'m',LineWidth=1.2)
        plot(t(2:end),y_plus(i,:),'g',LineWidth=1.2)
        plot(t(2:end),y_plus(i,:),'k',LineWidth=1.2)
    end
    ylabel(var{i},'Interpreter','latex')
    legend('Nominal', 'Perturbed','Noise','LKF','EKF')
end
xlabel('Time (secs)','Interpreter','latex')
sgtitle('Measurements vs Time','Interpreter','latex')

var = {'$e_{\xi_{g}}$ [m]','$e_{\eta_{g}}$ [m]','$e_{\theta_{g}}$ [rads]','$e_{\xi_{a}}$ [m]','$e_{\eta_{a}}$ [m]','$e_{\theta_{a}}$ [rads]'};
figure
for i = 1:n
    subplot(n,1,i); hold on;
    plot(t,error_x(i,:),'r',LineWidth=1.2)
    plot(t,2*sigma(i,:),'b--',LineWidth=1.2)
    plot(t,-2*sigma(i,:),'b--',LineWidth=1.2)
    if i == 3 || i == 6
        plot(t,wrapToPi(error_x(i,:)),'r',LineWidth=1.2)
        plot(t,wrapToPi(2*sigma(i,:)),'b--',LineWidth=1.2)
        plot(t,wrapToPi(-2*sigma(i,:)),'b--',LineWidth=1.2)
    else
        plot(t,error_x(i,:),'r',LineWidth=1.2)
        plot(t,2*sigma(i,:),'b--',LineWidth=1.2)
        plot(t,-2*sigma(i,:),'b--',LineWidth=1.2)
    end
    ylabel(var{i},'Interpreter','latex')
end
xlabel('Time (secs)','Interpreter','latex')
sgtitle('LKF States Error vs Time','Interpreter','latex')

var = {'$e_{\gamma_{ag}}$ [rads]','$e_{\rho_{g}}$ [m]','$e_{\gamma_{ga}}$ [rads]','$e_{\xi_{a}}$ [m]','$e_{\eta_{a}}$ [m]'};
figure
for i = 1:p
    subplot(p,1,i); hold on;
    if i == 1 || i == 3
        plot(t(2:end),wrapToPi(delta_y_minus(i,:)),'g',LineWidth=1.2)
        plot(t(2:end),wrapToPi(2*innov_cov(i,:)),'b--',LineWidth=1.2)
        plot(t(2:end),wrapToPi(-2*innov_cov(i,:)),'b--',LineWidth=1.2)
    else
        plot(t(2:end),delta_y_minus(i,:),'g',LineWidth=1.2)
        plot(t(2:end),2*innov_cov(i,:),'b--',LineWidth=1.2)
        plot(t(2:end),-2*innov_cov(i,:),'b--',LineWidth=1.2)
    end
    ylabel(var{i},'Interpreter','latex')
end
xlabel('Time (secs)','Interpreter','latex')
sgtitle('LKF Measurements Error vs Time','Interpreter','latex')

var = {'$e_{\xi_{g}}$ [m]','$e_{\eta_{g}}$ [m]','$e_{\theta_{g}}$ [rads]','$e_{\xi_{a}}$ [m]','$e_{\eta_{a}}$ [m]','$e_{\theta_{a}}$ [rads]'};
figure
for i = 1:n
    subplot(n,1,i); hold on;
    if i == 3 || i == 6
        plot(t,wrapToPi(error_xe(i,:)),'r',LineWidth=1.2)
        plot(t,wrapToPi(2*sigmae(i,:)),'b--',LineWidth=1.2)
        plot(t,wrapToPi(-2*sigmae(i,:)),'b--',LineWidth=1.2)
    else
        plot(t,error_xe(i,:),'r',LineWidth=1.2)
        plot(t,2*sigmae(i,:),'b--',LineWidth=1.2)
        plot(t,-2*sigmae(i,:),'b--',LineWidth=1.2)
    end
    ylabel(var{i},'Interpreter','latex')
end
xlabel('Time (secs)','Interpreter','latex')
sgtitle('EKF States Error vs Time','Interpreter','latex')

var = {'$e_{\gamma_{ag}}$ [rads]','$e_{\rho_{g}}$ [m]','$e_{\gamma_{ga}}$ [rads]','$e_{\xi_{a}}$ [m]','$e_{\eta_{a}}$ [m]'};
figure(21);
for i = 1:p
    subplot(p,1,i); hold on;
    if i == 1 || i == 3
        plot(t(2:end),wrapToPi(innovation(i,:)),'g',LineWidth=1.2)
        plot(t(2:end),wrapToPi(2*innov_cove(i,:)),'b--',LineWidth=1.2)
        plot(t(2:end),wrapToPi(-2*innov_cove(i,:)),'b--',LineWidth=1.2)
    else
        plot(t(2:end),innovation(i,:),'g',LineWidth=1.2)
        plot(t(2:end),2*innov_cove(i,:),'b--',LineWidth=1.2)
        plot(t(2:end),-2*innov_cove(i,:),'b--',LineWidth=1.2)
    end
    ylabel(var{i},'Interpreter','latex')
end
xlabel('Time (secs)','Interpreter','latex')
sgtitle('EKF Measurements Error vs Time','Interpreter','latex')

N = 10;

for i = 1:N
    [xnoise,ynoise] = TMTSim(t,x0,const,Qtrue,Rtrue,h);
    
    [x_plus, P_plus, sigma, y_plus, delta_y_minus, S, innov_cov] = LFK(delta_x0,P0,const,Q,Rtrue,xnom,ynom,ynoise);
    
    error_x = x_plus - xnoise;

    % testing NEES & NIS
    for ts = 1:length(t)
        invP(:,:,ts) = inv(P_plus(:,:,ts));
        error_x(3,ts) = wrapToPi(error_x(3,ts));
        error_x(6,ts) = wrapToPi(error_x(6,ts));
        NEES(i,ts) = error_x(:,ts)'*inv(P_plus(:,:,ts))*error_x(:,ts);
        if (ts >= 2)
            NIS(i,ts) = delta_y_minus(:,ts-1)'*inv(S(:,:,ts-1))*delta_y_minus(:,ts-1);
        end 
    end
end


% NEES Test:
epsNEESbar = mean(NEES,1);
alphaNEES = 0.05;
Nnx = N*n;
r1x = chi2inv(alphaNEES/2, Nnx)./N;
r2x = chi2inv(1-alphaNEES/2, Nnx)./N;

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
alphaNIS = 0.05;
Nny = N*p;
r1y = chi2inv(alphaNIS/2,Nny)./N;
r2y = chi2inv(1-alphaNIS/2,Nny)./N;

figure
plot(epsNISbar,'bo','MarkerSize',6,'LineWidth',2),hold on
plot(r1y*ones(size(epsNISbar)),'b--','LineWidth',2)
plot(r2y*ones(size(epsNISbar)),'b--','LineWidth',2)
ylabel('NIS statistic, $\bar{\epsilon}_y$','FontSize',14,'Interpreter','latex')
xlabel('time step, k','FontSize',14)
title('NIS Estimation Results','FontSize',14)
legend('NIS @ time k', 'r_1 bound', 'r_2 bound')
ylim([0 15])

N = 10;

for i = 1:N
    [xnoise,ynoise] = TMTSim(t,x0,const,Qtrue,Rtrue,h);
    
    [x_pluse, P_plus, sigmae, y_pluse, innovation, S, innov_cove] = EFK(x0,P0e,const,Qekf,Rtrue,ynoise,h);
    
    error_x = x_pluse - xnoise;

    % testing NEES & NIS
    for ts = 1:length(t)
        invP(:,:,ts) = inv(P_plus(:,:,ts));
        error_x(3,ts) = wrapToPi(error_x(3,ts));
        error_x(6,ts) = wrapToPi(error_x(6,ts));
        NEES(i,ts) = error_x(:,ts)'*inv(P_plus(:,:,ts))*error_x(:,ts);
        if (ts >= 2)
            NIS(i,ts) = innovation(:,ts-1)'*inv(S(:,:,ts-1))*innovation(:,ts-1);
        end 
    end
end


% NEES Test:
epsNEESbar = mean(NEES,1);
alphaNEES = 0.05;
Nnx = N*n;
r1x = chi2inv(alphaNEES/2, Nnx)./N;
r2x = chi2inv(1-alphaNEES/2, Nnx)./N;

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
alphaNIS = 0.01;
Nny = N*p;
r1y = chi2inv(alphaNIS/2,Nny)./N;
r2y = chi2inv(1-alphaNIS/2,Nny)./N;

figure
plot(epsNISbar,'bo','MarkerSize',6,'LineWidth',2),hold on
plot(r1y*ones(size(epsNISbar)),'b--','LineWidth',2)
plot(r2y*ones(size(epsNISbar)),'b--','LineWidth',2)
ylabel('NIS statistic, $\bar{\epsilon}_y$','FontSize',14,'Interpreter','latex')
xlabel('time step, k','FontSize',14)
title('NIS Estimation Results','FontSize',14)
legend('NIS @ time k', 'r_1 bound', 'r_2 bound')
ylim([0 15])

%% plot ekf
var = {'$\xi_{g}$ [m]','$\eta_{g}$ [m]','$\theta_{g}$ [rads]','$\xi_{a}$ [m]','$\eta_{a}$ [m]','$\theta_{a}$ [rads]'};
figure
for i = 1:n
    subplot(n,1,i)
    hold on
    if i == 3 || i == 6
        plot(t,wrapToPi(x_pluse(i,:)),'k',LineWidth=1.2)
        plot(t,wrapToPi(xnoise(i,:)),'m',LineWidth=1.2)
    else
        plot(t,x_pluse(i,:),'k',LineWidth=1.2)
        plot(t,xnoise(i,:),'m',LineWidth=1.2)
    end
    ylabel(var{i},'Interpreter','latex')
    legend('EKF', 'Noise')
end
xlabel('Time (secs)','Interpreter','latex')
sgtitle('EKF States vs Time with Noise','Interpreter','latex')

var = {'$\gamma_{ag}$ [rads]','$\rho_{g}$ [m]','$\gamma_{ga}$ [rads]','$\xi_{a}$ [m]','$\eta_{a}$ [m]'};
figure
for i = 1:p
    subplot(p,1,i)
    hold on
    if i == 1 || i == 3
        plot(t(2:end),wrapToPi(y_pluse(i,:)),'k',LineWidth=1.2)
        plot(t(2:end),wrapToPi(ynoise(i,:)),'m',LineWidth=1.2)
    else
        plot(t(2:end),y_pluse(i,:),'k',LineWidth=1.2)
        plot(t(2:end),ynoise(i,:),'m',LineWidth=1.2)
    end
    ylabel(var{i},'Interpreter','latex')
    legend('EKF', 'Noise')
end
xlabel('Time (secs)','Interpreter','latex')
sgtitle('EKF Measurements vs Time with Noise','Interpreter','latex')

%% Part 6
% Load necessary data and initialize
load('cooplocalization_finalproj_KFdata.mat');
[x_LKF, P_plus, sigma_LKF, y_plus, delta_y_minus, S, innov_cov] = LFK(delta_x0,P0,const,Q,Rtrue,xnom,ynom,ydata(:,2:end));
[x_EKF, P_plus, sigma_EKF, y_pluse, innovation, S, innov_cove] = EFK(x0,P0e,const,Qekf,Rtrue,ydata(:,2:end),h);

% Plot the results
figure;
for i = 1:6
    subplot(3, 2, i);
    hold on;
    if i == 3 || i == 6 
      plot(t, wrapToPi(x_LKF(i, :)), 'b', t, wrapToPi(x_EKF(i, :)), 'r');
      plot(t, wrapToPi(x_LKF(i, :)) + wrapToPi(2 * sigma_LKF(i, :)), 'b--', t, wrapToPi(x_LKF(i, :)) - wrapToPi(2 * sigma_LKF(i, :)), 'b--');
      plot(t, wrapToPi(x_EKF(i, :)) + wrapToPi(2 * sigma_EKF(i, :)), 'r--', t, wrapToPi(x_EKF(i, :) - 2 * sigma_EKF(i, :)), 'r--'); 
    else 
        plot(t, x_LKF(i, :), 'b', t, x_EKF(i, :), 'r');
        plot(t, x_LKF(i, :) + 2 * sigma_LKF(i, :), 'b--', t, x_LKF(i, :) - 2 * sigma_LKF(i, :), 'b--');
        plot(t, x_EKF(i, :) + 2 * sigma_EKF(i, :), 'r--', t, x_EKF(i, :) - 2 * sigma_EKF(i, :), 'r--');
    end 
    xlabel('Time [s]');
    ylabel(['State ', num2str(i)]);
    legend('LKF', 'EKF', 'LKF ±2σ', 'EKF ±2σ');
    title(['State ', num2str(i), ' Estimation' 'using Ydata']);
end

%%
function dx = dubinsEOM(t,x,w,const)
    theta_g = x(3);
    theta_a = x(6);
    xi_gdot = const.v_g*cos(theta_g) + w(1);
    eta_gdot = const.v_g*sin(theta_g) + w(2);
    theta_gdot = (const.v_g/const.L)*tan(const.phi_g) + w(3);
    xi_adot = const.v_a*cos(theta_a) + w(4);
    eta_adot = const.v_a*sin(theta_a) + w(5);
    theta_adot = const.omega_a + w(6);
    dx = [xi_gdot;eta_gdot;theta_gdot;xi_adot;eta_adot;theta_adot];
end

function y = calcY(x)
    % pull out all the states
    e_g(:,1) = x(:,1);
    n_g(:,1) = x(:,2);
    theta_g(:,1) = x(:,3);
    e_a(:,1) = x(:,4);
    n_a(:,1) = x(:,5);
    theta_a(:,1) = x(:,6);
    
    % calculate each row of the measurements
    y(1,:) = atan2((n_a-n_g),(e_a-e_g)) - theta_g;
    y(2,:) = sqrt((e_g-e_a).^2+(n_g-n_a).^2);
    y(3,:) = atan2((n_g-n_a),(e_g-e_a)) - theta_a;
    y(4,:) = e_a;
    y(5,:) = n_a;
end

function [x_plus, P_plus, sigma, y_plus, delta_y_minus, S, innov_cov] = LFK(delta_x0,P0,const,Q,R,xnom,ynom,ytrue)
    k = const.endTime/const.deltaT;
    delta_x_plus(:,1) = delta_x0;
    P_plus(:,:,1) = P0;
    sigma(:,1) = sqrt(diag(P_plus(:,:,1)));
    I = eye(length(delta_x0));
    for i = 1:k
        % Prediction Step
        [F_tilde,H_tilde] = findDTMatrices(xnom(i,:),const);
        delta_x_minus = F_tilde*delta_x_plus(:,i);
        P_minus = F_tilde*P_plus(:,:,i)*F_tilde' + Q;
        
        % Measurement Update
        S(:,:,i) = H_tilde*P_minus*H_tilde' + R;
        innov_cov(:,i) = sqrt(diag(S(:,:,i)));
        K = P_minus*H_tilde'*inv(S(:,:,i));
        P_plus(:,:,i+1) = (I-K*H_tilde)*P_minus;
        sigma(:,i+1) = sqrt(diag(P_plus(:,:,i+1)));
        delta_y_minus(:,i) = ytrue(:,i) - ynom(:,i);
        delta_y_minus(1,i) = wrapToPi(delta_y_minus(1,i));
        delta_y_minus(3,i) = wrapToPi(delta_y_minus(3,i));
        delta_x_plus(:,i+1) = delta_x_minus + K*(delta_y_minus(:,i)-H_tilde*delta_x_minus);
        delta_y_plus(:,i) = H_tilde*delta_x_plus(:,i+1);
    end
    x_plus = xnom' + delta_x_plus;
    y_plus = ynom + delta_y_plus;
end

function [x_plus, P_plus, sigma, y_plus, innovation, S, innov_cov] = EFK(x0,P0,const,Q,R,ytrue,h)
    k = const.endTime/const.deltaT;
    x_plus(:,1) = x0;
    P_plus(:,:,1) = P0;
    sigma(:,1) = sqrt(diag(P_plus(:,:,1)));
    I = eye(length(x0));
    for i = 1:k
        % Prediction Step
        time_dist = [0 const.deltaT];
        options = odeset('RelTol',1E-12,'AbsTol',1E-12);
        [t,X] = ode45(@(t,x) dubinsEOM(t,x,zeros(6,1),const),time_dist,x_plus(:,i),options); 
        x_minus = X(end,:)';
        [F_tilde,H_tilde] = findDTMatrices(x_minus,const);
        P_minus = F_tilde*P_plus(:,:,i)*F_tilde' + Q;
        
        % Measurement Update
        y_minus = h(x_minus);
        innovation(:,i) = ytrue(:,i) - y_minus;
        innovation(3,i) = wrapToPi(innovation(1,i));
        innovation(3,i) = wrapToPi(innovation(3,i));
        S(:,:,i) = H_tilde*P_minus*H_tilde' + R;
        innov_cov(:,i) = sqrt(diag(S(:,:,i)));
        K = P_minus*H_tilde'*inv(S(:,:,i));
        P_plus(:,:,i+1) = (I-K*H_tilde)*P_minus;
        sigma(:,i+1) = sqrt(diag(P_plus(:,:,i+1)));
        x_plus(:,i+1) = x_minus + K*innovation(:,i);
        y_plus(:,i) = h(x_plus(:,i+1));
    end
end

function [F,H] = findDTMatrices(x,const)
    xi_g = x(1);
    eta_g = x(2);
    theta_g = x(3);
    xi_a = x(4);
    eta_a = x(5);
    theta_a = x(6);
    F = eye(6) + 0.1*[0 0 -const.v_g*sin(theta_g) 0 0 0;
                      0 0  const.v_g*cos(theta_g) 0 0 0;
                      0 0  0                0 0 0;
                      0 0  0                0 0 -const.v_a*sin(theta_a);
                      0 0  0                0 0  const.v_a*cos(theta_a);
                      0 0  0                0 0 0];
  
    H = [ (eta_a - eta_g)/((xi_a - xi_g)^2 + (eta_a - eta_g)^2),  (xi_g - xi_a)/((xi_a - xi_g)^2 + (eta_a - eta_g)^2), -1,  (eta_g - eta_a)/((xi_a - xi_g)^2 + (eta_a - eta_g)^2),  (xi_a - xi_g)/((xi_a - xi_g)^2 + (eta_a - eta_g)^2),  0;
          (xi_g - xi_a)/sqrt((xi_g - xi_a)^2 + (eta_g - eta_a)^2),  (eta_g - eta_a)/sqrt((xi_g - xi_a)^2 + (eta_g - eta_a)^2),  0, (xi_a - xi_g)/sqrt((xi_g - xi_a)^2 + (eta_g - eta_a)^2),  (eta_a - eta_g)/sqrt((xi_g - xi_a)^2 + (eta_g - eta_a)^2),  0;
          (eta_a - eta_g)/((xi_g - xi_a)^2 + (eta_g - eta_a)^2),  (xi_g - xi_a)/((xi_g - xi_a)^2 + (eta_g - eta_a)^2),  0, (eta_g - eta_a)/((xi_g - xi_a)^2 + (eta_g - eta_a)^2),  (xi_a - xi_g)/((xi_g - xi_a)^2 + (eta_g - eta_a)^2), -1;
          0, 0, 0, 1, 0, 0;
          0, 0, 0, 0, 1, 0];
end

function [xnoise,ynoise] = TMTSim(t,x0,const,Qtrue,Rtrue,h)
    rng(100)
    xnoise(:,1) = x0;
    n = length(x0);
    p = 5;
    for i = 1: length(t)-1
        % Nonlinear Dynamics with perturbation + noise
        dt = [t(i) t(i+1)];
        options = odeset('RelTol',1E-12,'AbsTol',1E-12);
        w = mvnrnd(zeros(n,1),Qtrue)';
        [tnoise,Xnoise] = ode45(@(t,x) dubinsEOM(t,x,w,const),dt,xnoise(:,i),options);
        xnoise(:,i+1) = Xnoise(end,:)';
        v = mvnrnd(zeros(p,1),Rtrue)';
        ynoise(:,i) = h(xnoise(:,i+1)) + v;
    end
end

