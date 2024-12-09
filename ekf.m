% Initialize parameters
dt = 0.1; % Sampling time
N = 1000; % Number of time steps
MC_runs = 100; % Number of Monte Carlo simulations

% Initial state and covariance
x0 = [10; 0; pi/2; -60; 0; -pi/2]; % True initial state
x_hat = x0 + 0.1 * randn(6, 1); % Initial state estimate
P = eye(6); % Initial covariance

% Noise covariances
Q_EKF = diag([0.01, 0.01, 0.01, 0.01, 0.01, 0.01]);
R_EKF = diag([0.05, 0.1, 0.05, 0.01, 0.01]);

% Monte Carlo simulation
for mc = 1:MC_runs
    % Simulate true trajectory
    [x_true, y_meas] = simulate_trajectory(x0, Q_EKF, R_EKF, dt, N);
    
    % EKF loop
    for k = 1:N
        % Prediction step
        [x_hat_pred, F] = predict_state(x_hat, u_nom, dt);
        P_pred = F * P * F' + Q_EKF;
        
        % Update step
        [y_pred, H] = predict_measurement(x_hat_pred);
        S = H * P_pred * H' + R_EKF;
        K = P_pred * H' / S;
        x_hat = x_hat_pred + K * (y_meas(:, k) - y_pred);
        P = (eye(6) - K * H) * P_pred;
        
    end
end

