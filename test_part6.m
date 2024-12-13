% Load necessary data and initialize
load('cooplocalization_finalproj_KFdata.mat'); % Load provided observation data

% System constants and parameters
const.L = 0.5; % Length of the vehicle
const.v_g0 = 2; % Ground vehicle velocity
const.v_a0 = 12; % Aerial vehicle velocity
const.phi_g0 = -pi/18; % Steering angle
const.w_a0 = pi/25; % Angular rate
const.deltaT = 0.1; % Time step

% Initial conditions
x0 = zeros(6,1); % Initial state
P0 = eye(6); % Initial covariance matrix
Q = 0.01 * eye(6); % Process noise covariance
R = Rtrue; % Use provided measurement noise covariance

% Nominal trajectory
endTime = tvec(end); % Use time vector from provided data
t = tvec;
x_nom = zeros(length(t), 6);
for i = 1:length(t)
    x_nom(i,:) = FindNominal(t(i), x0, const.v_g0, const.v_a0, const.L, const.phi_g0, const.w_a0);
end

y_nom = findYnom(x_nom); % Calculate nominal measurements

% Observation data
%y_meas = ydata'; % Provided observation data

% Run Linearized Kalman Filter (LKF)
%[x_LKF, P_LKF, ~, y_LKF, sigma_LKF, ~] = LKF(x0, P0, const, @CT_to_DT, x_nom, y_nom, ydata, Q, R); % changed y_meas to ydata
[x_plus, P_plus, sigma, y_plus, delta_y_minus, S, innov_cov] = LFK(delta_x0,P0,const,Q,Rtrue,xnom,ynom,ydata);

% Run Extended Kalman Filter (EKF)
[x_EKF, P_EKF, ~, ~, ~, ~, sigma_EKF] = ekf(y_meas, x0, P0, Q, R, zeros(4, length(t)), const.deltaT, length(t), @FindNominal, @findYnom, @CT_to_DT, @findYnom);

% Plot the results
figure;
for i = 1:6
    subplot(3, 2, i);
    plot(t, x_LKF(i, :), 'b', t, x_EKF(i, :), 'r');
    hold on;
    plot(t, x_LKF(i, :) + 2 * sigma_LKF(i, :), 'b--', t, x_LKF(i, :) - 2 * sigma_LKF(i, :), 'b--');
    plot(t, x_EKF(i, :) + 2 * sigma_EKF(i, :), 'r--', t, x_EKF(i, :) - 2 * sigma_EKF(i, :), 'r--');
    xlabel('Time [s]');
    ylabel(['State ', num2str(i)]);
    legend('LKF', 'EKF', 'LKF ±2σ', 'EKF ±2σ');
    title(['State ', num2str(i), ' Estimation']);
end

% Compare performance (e.g., RMSE or NEES/NIS metrics)
[~, ~, ~, ~, ~, NEES_LKF, NIS_LKF] = FindNISNESS(1, x0, P0, x_nom, y_nom, @CT_to_DT, const, Q, R, endTime);
[~, ~, ~, ~, ~, NEES_EKF, NIS_EKF] = FindNISNESS(1, x0, P0, x_nom, y_nom, @CT_to_DT, const, Q, R, endTime);

fprintf('Average NEES (LKF): %.2f\n', mean(NEES_LKF));
fprintf('Average NEES (EKF): %.2f\n', mean(NEES_EKF));
fprintf('Average NIS (LKF): %.2f\n', mean(NIS_LKF));
fprintf('Average NIS (EKF): %.2f\n', mean(NIS_EKF));
