% Initialization
dt = 0.1; % Sampling time
N = 1000; % Number of time steps
MC_runs = 100; % Number of Monte Carlo simulations

% Load initial parameters from prog_report_1
x0 = [10; 0; pi/2; -60; 0; -pi/2]; % Initial state
u_nom = [2; -pi/18; 12; pi/25]; % Control inputs

% Process and measurement noise covariances
Q_EKF = diag([0.01, 0.01, 0.01, 0.01, 0.01, 0.01]);
R_EKF = diag([0.05, 0.1, 0.05, 0.01, 0.01]);


% Monte Carlo Simulations
for mc = 1:MC_runs
    % Generate nominal trajectory and noisy ground truth
    [t, x_nom] = FindNominal(x0, u_nom, dt, N); % Nominal trajectory
    w = mvnrnd(zeros(6, 1), Q_EKF, N)'; % Process noise
    x_true = x_nom + w; % Noisy ground truth trajectory
    
    % Generate noisy measurements
    y_nom = findYnom(x_nom); % Nominal measurements
    v = mvnrnd(zeros(5, 1), R_EKF, N)'; % Measurement noise
    y_meas = y_nom + v; % Noisy measurements
    
    % EKF Initialization
    x_hat = x0; % Initial state estimate
    P = eye(6); % Initial covariance
    
    % EKF Loop
    for k = 1:N
        % Prediction Step
        [x_hat_pred, F] = ugvEOM(x_hat, u_nom, dt); % Predict state and Jacobian
        P_pred = F * P * F' + Q_EKF; % Predict covariance
        
        % Measurement Update
        [y_pred, H] = findYnom(x_hat_pred); % Predicted measurement and Jacobian
        S = H * P_pred * H' + R_EKF; % Innovation covariance
        K = P_pred * H' / S; % Kalman gain
        x_hat = x_hat_pred + K * (y_meas(:, k) - y_pred); % Update state
        P = (eye(6) - K * H) * P_pred; % Update covariance
        
    end
end

