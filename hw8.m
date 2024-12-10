%5044 Hw8

clear; 
clc;
close all;

% Q1
del_t = 0.5; 
omega_A = 0.045; %rad
omega_B = -0.045;

F =  @(omega) [1 (sin(omega*del_t)/omega) 0 -(1-cos(omega*del_t))/omega;
    0 cos(omega*del_t) 0 -sin(omega*del_t);
    0 (1-cos(omega*del_t))/omega 1 (sin(omega*del_t)/omega);
    0 sin(omega*del_t) 0 cos(omega*del_t)];

n = 4; %num states
Fa = F(omega_A);
Fb = F(omega_B);

Gamma = [0 0; 1 0; 0 0; 0 1];
qw = 10; %m/s^2
W = qw*[2 0.05; 0.05 0.5];
% finding Q
A = @(omega) [0 1 0 0; 0 0 0 -omega; 0 0 0 1; 0 omega 0 0];
gwg = Gamma*W*Gamma';
Za = del_t * [[-A(omega_A);zeros(n)], [gwg; A(omega_A)']];
Zb = del_t * [[-A(omega_B);zeros(n)], [gwg; A(omega_B)']];
% matrix exponential of the Z
ez_a = expm(Za);
ez_b = expm(Zb);
% pull out Q from exp Z
%Fa = ez_a(n+1:end,n+1:end)'
Qa = Fa * ez_a(1:n,n+1:end);
Qb = Fb * ez_b(1:n,n+1:end);


%% Q2

% a 
% fix random seed
rng(100);

load('hw8problemdata.mat')

Ra = [20 0.05; 0.05 20];
%Ra = [8000 500; 500 8000];
H = [1 0 0 0; 0 0 1 0];

T = length(xasingle_truth); %number of samples
p = length(Ra); %num measurements in each sample
Sv = chol(Ra,'lower');

%create vector of gaussian random variables
qk = mvnrnd(zeros(p,1),eye(p),T);
%qk = randn(length(xasingle_truth),1)

% create noisy measurements Hx + v (v = Sv*qk)
y_noisy = (H * xasingle_truth + Sv * qk')';

k_time = 1:T; 

% plot components of ya vs time for 1st 20 seconds
figure(); hold on;
subplot(2,1,1); hold on;
plot(k_time(1:40)*del_t,y_noisy(1:40,1),LineWidth=1.3)
ylabel('Measured East_A')
subplot(2,1,2); hold on;
plot(k_time(1:40)*del_t,y_noisy(1:40,2),LineWidth=1.3)
ylabel('Measured North_A')
xlabel('Time (s)')
sgtitle('Simulated Noisy Measurements y_A vs Time')

% b 
ua0 = [0; 85*cos(pi/4); 0; -85*sin(pi/4)];
Pa0 = 900*diag([10, 2, 10, 2]);

%%% Kalman 
% initialize
x_plus(:,1) = ua0;
P_plus(:,:,1) = Pa0;

I = eye(n);
H_kal = repmat(H,[1,1,T+1]);
Q = Qa;
R = repmat(Ra,[1,1,T+1]);

for k = 1:T

    %%%%%%%%%%%%%%%%%%%%%%%%
    %%% prediction step section 
    x_minus(:,k+1) = Fa * x_plus(:,k);
    P_minus(:,:,k+1) = Fa * P_plus(:,:,k) * Fa' + Q;
    K(:,:,k+1) = P_minus(:,:,k+1) * H_kal(:,:,k+1)' * inv(H_kal(:,:,k+1)*P_minus(:,:,k+1)*H_kal(:,:,k+1)'+R(:,:,k+1));

    % % transition from minus to plus (testing the prediciton step)
    % x_plus(:,k+1) = x_minus(:,k+1);
    % P_plus(:,:,k+1) = P_minus(:,:,k+1);

    %%%%%%%%%%%%%%%%%%%%%%%
    %%% correction step 
    x_plus(:,k+1) = x_minus(:,k+1) + K(:,:,k+1) * (y_noisy(k,:)' - H_kal(:,:,k+1) * x_minus(:,k+1));
    P_plus(:,:,k+1) = (I - K(:,:,k+1) * H_kal(:,:,k+1)) * P_minus(:,:,k+1);

end

% plotting check 
figure(); 
subplot(4,1,1); hold on;
plot([0,k_time],x_plus(1,:))
subplot(4,1,2); hold on;
plot([0,k_time],x_plus(2,:))
subplot(4,1,3); hold on;
plot([0,k_time],x_plus(3,:))
subplot(4,1,4); hold on;
plot([0,k_time],x_plus(4,:))
title('States vs Time')


% plot estimated state error vs time
figure(); 
subplot(4,1,1); hold on;
plot(k_time,x_plus(1,2:end)-xasingle_truth(1,:))
ylabel('$\xi_a$ [m]', 'Interpreter','latex')
xlim([0 202])
subplot(4,1,2); hold on;
plot(k_time,x_plus(2,2:end)-xasingle_truth(2,:))
ylabel('$\dot{\xi_a}$ [m/s]', 'Interpreter','latex')
xlim([0 202])
subplot(4,1,3); hold on;
plot(k_time,x_plus(3,2:end)-xasingle_truth(3,:))
ylabel('$\eta_a$ [m]', 'Interpreter','latex')
xlim([0 202])
subplot(4,1,4); hold on;
plot(k_time,x_plus(4,2:end)-xasingle_truth(4,:))
ylabel('$\dot{\eta_a}$ [m/s]', 'Interpreter','latex')
xlim([0 202])
xlabel('time (k)')
sgtitle('State Errors vs Time (2b)')

% plot estimated state error vs time with sigma bounds
figure(); 
subplot(4,1,1); hold on;
plot(k_time,x_plus(1,2:end)-xasingle_truth(1,:),'k',LineWidth=1.5)
plot(k_time,+ 2*sqrt(reshape(P_plus(1,1,2:end),[1,length(P_plus)-1])),'b--',LineWidth=1.2)
plot(k_time,- 2*sqrt(reshape(P_plus(1,1,2:end),[1,length(P_plus)-1])),'b--',LineWidth=1.2)
ylabel('$\xi_a error$ [m]', 'Interpreter','latex')
xlim([0 202])
subplot(4,1,2); hold on;
plot(k_time,x_plus(2,2:end)-xasingle_truth(2,:),'k',LineWidth=1.5)
plot(k_time,+ 2*sqrt(reshape(P_plus(2,2,2:end),[1,length(P_plus)-1])),'b--',LineWidth=1.2)
plot(k_time,- 2*sqrt(reshape(P_plus(2,2,2:end),[1,length(P_plus)-1])),'b--',LineWidth=1.2)
ylabel('$\dot{\xi_a} error$ [m/s]', 'Interpreter','latex')
xlim([0 202])
subplot(4,1,3); hold on;
plot(k_time,x_plus(3,2:end)-xasingle_truth(3,:),'k',LineWidth=1.5)
plot(k_time,+ 2*sqrt(reshape(P_plus(3,3,2:end),[1,length(P_plus)-1])),'b--',LineWidth=1.2)
plot(k_time,- 2*sqrt(reshape(P_plus(3,3,2:end),[1,length(P_plus)-1])),'b--',LineWidth=1.2)
ylabel('$\eta_a error$ [m]', 'Interpreter','latex')
xlim([0 202])
subplot(4,1,4); hold on;
plot(k_time,x_plus(4,2:end)-xasingle_truth(4,:),'k',LineWidth=1.5)
plot(k_time,+ 2*sqrt(reshape(P_plus(4,4,2:end),[1,length(P_plus)-1])),'b--',LineWidth=1.2)
plot(k_time,- 2*sqrt(reshape(P_plus(4,4,2:end),[1,length(P_plus)-1])),'b--',LineWidth=1.2)
ylabel('$\dot{\eta_a} error$ [m/s]', 'Interpreter','latex')
xlim([0 202])
xlabel('time (k)')
sgtitle('State Errors vs Time with 2 Sigma Bounds (2b)')

% plot estimated 2 sigma bounds vs time
figure(); 
subplot(4,1,1); hold on;
plot([0,k_time],2*sqrt(reshape(P_plus(1,1,:),[1,length(P_plus)])))
ylabel('$2\sigma \xi_a$ [m]', 'Interpreter','latex')
xlim([0 202])
subplot(4,1,2); hold on;
plot([0,k_time],2*sqrt(reshape(P_plus(2,2,:),[1,length(P_plus)])))
ylabel('$2\sigma \dot{\xi_a}$ [m/s]', 'Interpreter','latex')
xlim([0 202])
subplot(4,1,3); hold on;
plot([0,k_time],2*sqrt(reshape(P_plus(3,3,:),[1,length(P_plus)])))
ylabel('$2\sigma \eta_a$ [m]', 'Interpreter','latex')
xlim([0 202])
subplot(4,1,4); hold on;
plot([0,k_time],2*sqrt(reshape(P_plus(4,4,:),[1,length(P_plus)])))
ylabel('$2\sigma \dot{\eta_a}$ [m/s]', 'Interpreter','latex')
xlim([0 202])
xlabel('time (k)')
sgtitle('2 Sigma Bounds vs Time')

% checking plots by plotting states and meas together 
figure(); hold on;
subplot(2,1,1); hold on;
plot([0,k_time],x_plus(1,:),'b')
plot(k_time,y_noisy(:,1),'r--')
ylabel('Easting A')
subplot(2,1,2); hold on;
plot([0,k_time],x_plus(3,:),'b')
plot(k_time,y_noisy(:,2),'r--')
ylabel('Northing A')
xlabel('Time (k)')
sgtitle('checking 2')

%% Q3
% initial values 
ua0 = [0; 85*cos(pi/4); 0; -85*sin(pi/4)];
Pa0 = 900*diag([10, 2, 10, 2]);
ub0 = [3200; 85*cos(pi/4); 3200; -85*sin(pi/4)];
Pb0 = 900*diag([11, 4, 11, 4]);

Ra = [20 0.05; 0.05 20];
Rd = [10 0.15; 0.15 10];
Ha = [1 0 0 0; 0 0 1 0];
Hd = [1 0 0 0 -1 0 0 0;
      0 0 1 0 0 0 -1 0];

% simulate ya and yd noisy meas
T = length(xadouble_truth); %number of samples

% get ya noisy
pa = length(Ra); %num measurements in each sample
Sva = chol(Ra,'lower');
%create vector of gaussian random variables
qka = mvnrnd(zeros(pa,1),eye(pa),T);
% create noisy measurements Hx + v (v = Sv*qk)
ya_noisy = (Ha * xadouble_truth + Sva * qka')';

% get yd noisy
pd = length(Rd); %num measurements in each sample
Svd = chol(Rd,'lower');
%create vector of gaussian random variables
qkd = mvnrnd(zeros(pd,1),eye(pd),T);
% create noisy measurements Hx + v (v = Sv*qk)
yd_noisy = (Hd * [xadouble_truth; xbdouble_truth] + Svd * qkd')';

%%% a -- use ya and yd 
% making all the combined matricies
F_3 = blkdiag(Fa, Fb);
R_3 = blkdiag(Ra, Rd);
P0_3 = blkdiag(Pa0, Pb0);
Q_3 = blkdiag(Qa, Qb);
H_3 = [1 0 0 0 0 0 0 0;
       0 0 1 0 0 0 0 0;
       1 0 0 0 -1 0 0 0;
       0 0 1 0 0 0 -1 0];
u0_3 = [ua0; ub0];

% run the kalman filter with my new function 
[x_plus_3, P_plus_3, innovation_3] = Kalman_Filter(u0_3, P0_3, F_3, H_3, R_3, Q_3, [ya_noisy';yd_noisy']);

%%% plot results 
time_k = 1:length(ya_noisy);
% 
%   plot states with measurements vs time 
    % figure(); hold on;
    % subplot(4,1,1); hold on;
    % plot([0,time_k],x_plus_3(1,:),'b')
    % plot(time_k,ya_noisy(:,1),'r--')
    % ylabel('East A')
    % subplot(4,1,2); hold on;
    % plot([0,time_k],x_plus_3(3,:),'b')
    % plot(time_k,ya_noisy(:,2),'r--')
    % ylabel('North A')
    % subplot(4,1,3); hold on;
    % plot([0,time_k],x_plus_3(1,:)-x_plus_3(5,:),'b')
    % plot(time_k,yd_noisy(:,1),'r--')
    % ylabel('East A - East B')
    % subplot(4,1,4); hold on;
    % plot([0,time_k],x_plus_3(3,:)-x_plus_3(7,:),'b')
    % plot(time_k,yd_noisy(:,2),'r--')
    % ylabel('North A - North B')
    % xlabel('Time (k)')
    % sgtitle('States and Measurements vs Time')
    % 
    % % plot the covaraiances 
    % figure(); hold on;
    % subplot(4,2,1); hold on;
    % plot([0,time_k],2*sqrt(reshape(P_plus_3(1,1,:),[1,length(P_plus_3)])),'b')
    % ylabel('East A Cov')
    % subplot(4,2,2); hold on;
    % plot([0,time_k],2*sqrt(reshape(P_plus_3(2,2,:),[1,length(P_plus_3)])),'b')
    % ylabel('East_dot A Cov')
    % subplot(4,2,3); hold on;
    % plot([0,time_k],2*sqrt(reshape(P_plus_3(3,3,:),[1,length(P_plus_3)])),'b')
    % ylabel('North A Cov')
    % subplot(4,2,4); hold on;
    % plot([0,time_k],2*sqrt(reshape(P_plus_3(5,5,:),[1,length(P_plus_3)])),'b')
    % ylabel('North_dot A Cov')
    % subplot(4,2,5); hold on;
    % plot([0,time_k],2*sqrt(reshape(P_plus_3(6,6,:),[1,length(P_plus_3)])),'b')
    % ylabel('East B Cov')
    % subplot(4,2,6); hold on;
    % plot([0,time_k],2*sqrt(reshape(P_plus_3(7,7,:),[1,length(P_plus_3)])),'b')
    % ylabel('East_dot B Cov')
    % subplot(4,2,7); hold on;
    % plot([0,time_k],2*sqrt(reshape(P_plus_3(8,8,:),[1,length(P_plus_3)])),'b')
    % ylabel('North B Cov')
    % subplot(4,2,8); hold on;
    % plot([0,time_k],2*sqrt(reshape(P_plus_3(4,4,:),[1,length(P_plus_3)])),'b')
    % ylabel('North_dot B Cov')
    % xlabel('Time (k)')
    % sgtitle('Covariances vs Time')
    % 
    % % plot position errors vs time 
    % % plot estimated state error vs time
    % figure(); 
    % subplot(4,1,1); hold on;
    % plot(k_time,x_plus_3(1,2:end)-xadouble_truth(1,:))
    % ylabel('$\xi_a$ [m]', 'Interpreter','latex')
    % xlim([0 202])
    % subplot(4,1,2); hold on;
    % plot(k_time,x_plus_3(2,2:end)-xadouble_truth(2,:))
    % ylabel('$\dot{\xi_a}$ [m/s]', 'Interpreter','latex')
    % xlim([0 202])
    % subplot(4,1,3); hold on;
    % plot(k_time,x_plus_3(3,2:end)-xadouble_truth(3,:))
    % ylabel('$\eta_a$ [m]', 'Interpreter','latex')
    % xlim([0 202])
    % subplot(4,1,4); hold on;
    % plot(k_time,x_plus_3(4,2:end)-xadouble_truth(4,:))
    % ylabel('$\dot{\eta_a}$ [m/s]', 'Interpreter','latex')
    % xlim([0 202])
    % xlabel('time (k)')
    % sgtitle('State A Errors vs Time')
    % 
    % figure(); 
    % subplot(4,1,1); hold on;
    % plot(k_time,x_plus_3(5,2:end)-xbdouble_truth(1,:))
    % ylabel('$\xi_a$ [m]', 'Interpreter','latex')
    % xlim([0 202])
    % subplot(4,1,2); hold on;
    % plot(k_time,x_plus_3(6,2:end)-xbdouble_truth(2,:))
    % ylabel('$\dot{\xi_a}$ [m/s]', 'Interpreter','latex')
    % xlim([0 202])
    % subplot(4,1,3); hold on;
    % plot(k_time,x_plus_3(7,2:end)-xbdouble_truth(3,:))
    % ylabel('$\eta_a$ [m]', 'Interpreter','latex')
    % xlim([0 202])
    % subplot(4,1,4); hold on;
    % plot(k_time,x_plus_3(8,2:end)-xbdouble_truth(4,:))
    % ylabel('$\dot{\eta_a}$ [m/s]', 'Interpreter','latex')
    % xlim([0 202])
    % xlabel('time (k)')
    % sgtitle('State B Errors vs Time')
%

% plot estimated state error vs time with sigma bounds
figure(); 
subplot(2,1,1); hold on;
plot(k_time,x_plus_3(1,2:end)-xadouble_truth(1,:),'k',LineWidth=1.5)
plot(k_time,+ 2*sqrt(reshape(P_plus_3(1,1,2:end),[1,length(P_plus_3)-1])),'b--',LineWidth=1.2)
plot(k_time,- 2*sqrt(reshape(P_plus_3(1,1,2:end),[1,length(P_plus_3)-1])),'b--',LineWidth=1.2)
ylabel('$\xi_a error$ [m]', 'Interpreter','latex')
xlim([0 202])
subplot(2,1,2); hold on;
plot(k_time,x_plus_3(3,2:end)-xadouble_truth(3,:),'k',LineWidth=1.5)
plot(k_time,+ 2*sqrt(reshape(P_plus_3(3,3,2:end),[1,length(P_plus_3)-1])),'b--',LineWidth=1.2)
plot(k_time,- 2*sqrt(reshape(P_plus_3(3,3,2:end),[1,length(P_plus_3)-1])),'b--',LineWidth=1.2)
ylabel('$\eta_a error$ [m]', 'Interpreter','latex')
xlim([0 202])
xlabel('time (k)')
sgtitle('State A Errors vs Time with 2 Sigma Bounds (3a)')


figure(); 
subplot(2,1,1); hold on;
plot(k_time,x_plus_3(5,2:end)-xbdouble_truth(1,:),'k',LineWidth=1.5)
plot(k_time,+ 2*sqrt(reshape(P_plus_3(5,5,2:end),[1,length(P_plus_3)-1])),'b--',LineWidth=1.2)
plot(k_time,- 2*sqrt(reshape(P_plus_3(5,5,2:end),[1,length(P_plus_3)-1])),'b--',LineWidth=1.2)
ylabel('$\xi_b error$ [m]', 'Interpreter','latex')
xlim([0 202])
subplot(2,1,2); hold on;
plot(k_time,x_plus_3(7,2:end)-xbdouble_truth(3,:),'k',LineWidth=1.5)
plot(k_time,+ 2*sqrt(reshape(P_plus_3(7,7,2:end),[1,length(P_plus_3)-1])),'b--',LineWidth=1.2)
plot(k_time,- 2*sqrt(reshape(P_plus_3(7,7,2:end),[1,length(P_plus_3)-1])),'b--',LineWidth=1.2)
ylabel('$\eta_b error$ [m]', 'Interpreter','latex')
xlim([0 202])
xlabel('time (k)')
sgtitle('State B Errors vs Time with 2 Sigma Bounds (3a)')

%%% b -- use yd only 
ua0 = [0; 85*cos(pi/4); 0; -85*sin(pi/4)];
Pa0 = 900*diag([10, 2, 10, 2]);
ub0 = [3200; 85*cos(pi/4); 3200; -85*sin(pi/4)];
Pb0 = 900*diag([11, 4, 11, 4]);
Rd = [10 0.15; 0.15 10];
Hd = [1 0 0 0 -1 0 0 0;
      0 0 1 0 0 0 -1 0];
F_3 = blkdiag(Fa, Fb);
P0_3 = blkdiag(Pa0, Pb0);
Q_3 = blkdiag(Qa, Qb);
u0_3 = [ua0; ub0];

% run the kalman filter with my new function 
[x_plus_3b, P_plus_3b, innovation_3b] = Kalman_Filter(u0_3, P0_3, F_3, Hd, Rd, Q_3, yd_noisy');

%%% plot results 
time_k = 1:length(yd_noisy);
% %plot states with measurements vs time 
% figure(); hold on;
% subplot(2,1,1); hold on;
% plot([0,time_k],x_plus_3b(1,:)-x_plus_3b(5,:),'b')
% plot(time_k,yd_noisy(:,1),'r--')
% ylabel('East A - East B')
% subplot(2,1,2); hold on;
% plot([0,time_k],x_plus_3b(3,:)-x_plus_3b(7,:),'b')
% plot(time_k,yd_noisy(:,2),'r--')
% ylabel('North A - North B')
% xlabel('Time (k)')
% sgtitle('States and Measurements vs Time')
% 
% % plot estimated state error vs time
% figure(); 
% subplot(4,1,1); hold on;
% plot(k_time,x_plus_3b(1,2:end)-xadouble_truth(1,:))
% ylabel('$\xi_a$ [m]', 'Interpreter','latex')
% xlim([0 202])
% subplot(4,1,2); hold on;
% plot(k_time,x_plus_3b(2,2:end)-xadouble_truth(2,:))
% ylabel('$\dot{\xi_a}$ [m/s]', 'Interpreter','latex')
% xlim([0 202])
% subplot(4,1,3); hold on;
% plot(k_time,x_plus_3b(3,2:end)-xadouble_truth(3,:))
% ylabel('$\eta_a$ [m]', 'Interpreter','latex')
% xlim([0 202])
% subplot(4,1,4); hold on;
% plot(k_time,x_plus_3b(4,2:end)-xadouble_truth(4,:))
% ylabel('$\dot{\eta_a}$ [m/s]', 'Interpreter','latex')
% xlim([0 202])
% xlabel('time (k)')
% sgtitle('State A Errors vs Time')
% 
% figure(); 
% subplot(4,1,1); hold on;
% plot(k_time,x_plus_3b(5,2:end)-xbdouble_truth(1,:))
% ylabel('$\xi_a$ [m]', 'Interpreter','latex')
% xlim([0 202])
% subplot(4,1,2); hold on;
% plot(k_time,x_plus_3b(6,2:end)-xbdouble_truth(2,:))
% ylabel('$\dot{\xi_a}$ [m/s]', 'Interpreter','latex')
% xlim([0 202])
% subplot(4,1,3); hold on;
% plot(k_time,x_plus_3b(7,2:end)-xbdouble_truth(3,:))
% ylabel('$\eta_a$ [m]', 'Interpreter','latex')
% xlim([0 202])
% subplot(4,1,4); hold on;
% plot(k_time,x_plus_3b(8,2:end)-xbdouble_truth(4,:))
% ylabel('$\dot{\eta_a}$ [m/s]', 'Interpreter','latex')
% xlim([0 202])
% xlabel('time (k)')
% sgtitle('State B Errors vs Time')

% plot estimated state error vs time with sigma bounds
figure(); 
subplot(2,1,1); hold on;
plot(k_time,x_plus_3b(1,2:end)-xadouble_truth(1,:),'k',LineWidth=1.5)
plot(k_time,+ 2*sqrt(reshape(P_plus_3b(1,1,2:end),[1,length(P_plus_3b)-1])),'b--',LineWidth=1.2)
plot(k_time,- 2*sqrt(reshape(P_plus_3b(1,1,2:end),[1,length(P_plus_3b)-1])),'b--',LineWidth=1.2)
ylabel('$\xi_a error$ [m]', 'Interpreter','latex')
xlim([0 202])
subplot(2,1,2); hold on;
plot(k_time,x_plus_3b(3,2:end)-xadouble_truth(3,:),'k',LineWidth=1.5)
plot(k_time,+ 2*sqrt(reshape(P_plus_3b(3,3,2:end),[1,length(P_plus_3b)-1])),'b--',LineWidth=1.2)
plot(k_time,- 2*sqrt(reshape(P_plus_3b(3,3,2:end),[1,length(P_plus_3b)-1])),'b--',LineWidth=1.2)
ylabel('$\eta_a error$ [m]', 'Interpreter','latex')
xlim([0 202])
xlabel('time (k)')
sgtitle('State A Errors vs Time with 2 Sigma Bounds (3b)')

figure(); 
subplot(2,1,1); hold on;
plot(k_time,x_plus_3b(5,2:end)-xbdouble_truth(1,:),'k',LineWidth=1.5)
plot(k_time,+ 2*sqrt(reshape(P_plus_3b(5,5,2:end),[1,length(P_plus_3b)-1])),'b--',LineWidth=1.2)
plot(k_time,- 2*sqrt(reshape(P_plus_3b(5,5,2:end),[1,length(P_plus_3b)-1])),'b--',LineWidth=1.2)
ylabel('$\xi_b error$ [m]', 'Interpreter','latex')
xlim([0 202])
subplot(2,1,2); hold on;
plot(k_time,x_plus_3b(7,2:end)-xbdouble_truth(3,:),'k',LineWidth=1.5)
plot(k_time,+ 2*sqrt(reshape(P_plus_3b(7,7,2:end),[1,length(P_plus_3b)-1])),'b--',LineWidth=1.2)
plot(k_time,- 2*sqrt(reshape(P_plus_3b(7,7,2:end),[1,length(P_plus_3b)-1])),'b--',LineWidth=1.2)
ylabel('$\eta_b error$ [m]', 'Interpreter','latex')
xlim([0 202])
xlabel('time (k)')
sgtitle('State B Errors vs Time with 2 Sigma Bounds (3b)')



%% Testing NEES and NIS
n = 8;
p = 4;
N = 1000;
for i = 1:N
    % get ya noisy
    pa = length(Ra); %num measurements in each sample
    Sva = chol(Ra,'lower');
    %create vector of gaussian random variables
    qka = mvnrnd(zeros(pa,1),eye(pa),T);
    % create noisy measurements Hx + v (v = Sv*qk)
    ya_noisy = (Ha * xadouble_truth + Sva * qka')';
    
    % get yd noisy
    pd = length(Rd); %num measurements in each sample
    Svd = chol(Rd,'lower');
    %create vector of gaussian random variables
    qkd = mvnrnd(zeros(pd,1),eye(pd),T);
    % create noisy measurements Hx + v (v = Sv*qk)
    yd_noisy = (Hd * [xadouble_truth; xbdouble_truth] + Svd * qkd')';
    
    y_truth = [ya_noisy; yd_noisy];

    %%% a -- use ya and yd 
    % making all the combined matricies
    F_3 = blkdiag(Fa, Fb);
    R_3 = blkdiag(Ra, Rd);
    P0_3 = blkdiag(Pa0, Pb0);
    Q_3 = blkdiag(Qa, Qb);
    H_3 = [1 0 0 0 0 0 0 0;
           0 0 1 0 0 0 0 0;
           1 0 0 0 -1 0 0 0;
           0 0 1 0 0 0 -1 0];
    u0_3 = [ua0; ub0];
    
    % run the kalman filter with my new function 
    [x_plus_3, P_plus_3, innovation_3] = Kalman_Filter(u0_3, P0_3, F_3, H_3, R_3, Q_3, [ya_noisy';yd_noisy']);


    % Calculate NEES and NIS
    for j = 1: length(t)
        NEES(i,j) = (xtrue(j,:) - x_plus(:,j)')*inv(P_plus(:,:,j))*(xtrue(j,:) - x_plus(:,j)')';
        if j > 2
            NIS(i,j-1) = (ytrue(:,j-1) - y_minus(:,j-1))'*inv(Sk(:,:,j-1))*(ytrue(:,j-1) - y_minus(:,j-1));
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
ylabel('NEES statistic, \bar{\epsilon}_x','FontSize',14)
xlabel('time step, k','FontSize',14)
title('NEES Estimation Results','FontSize',14)
legend('NEES @ time k', 'r_1 bound', 'r_2 bound')

epsNISbar = mean(NIS,1);
alphaNIS = 0.05;
Nny = N*p;
r1y = chi2inv(alphaNIS/2,Nny)./N;
r2y = chi2inv(1-alphaNIS/2,Nny)./N;

figure
plot(epsNISbar,'bo','MarkerSize',6,'LineWidth',2),hold on
plot(r1y*ones(size(epsNISbar)),'b--','LineWidth',2)
plot(r2y*ones(size(epsNISbar)),'b--','LineWidth',2)
ylabel('NIS statistic, \bar{\epsilon}_y','FontSize',14)
xlabel('time step, k','FontSize',14)
title('NIS Estimation Results','FontSize',14)
legend('NIS @ time k', 'r_1 bound', 'r_2 bound')


