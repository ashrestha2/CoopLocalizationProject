function [x_plus, P_plus, innovations, F_matrices] = ekf(y_meas, x0, P0, Q, R, u, dt, N, f, h, F_func, H_func)
    % Extended Kalman Filter (EKF) with F_k output
    %
    % Outputs:
    % x_plus      - Updated state estimates (state_dim x N)
    % P_plus      - Updated covariance matrices (state_dim x state_dim x N)
    % innovations - Measurement residuals (meas_dim x N)
    % F_matrices  - State transition Jacobians (state_dim x state_dim x N)

    % Initialize dimensions
    state_dim = length(x0);
    meas_dim = size(y_meas, 1);

    % Preallocate storage
    x_plus = zeros(state_dim, N); % Posterior state estimates
    P_plus = zeros(state_dim, state_dim, N); % Posterior covariance matrices
    innovations = zeros(meas_dim, N); % Measurement residuals (innovations)
    F_matrices = zeros(state_dim, state_dim, N); % State transition Jacobians

    % Initial conditions
    x_plus(:, 1) = x0; % Initial state estimate
    P_plus(:, :, 1) = P0; % Initial covariance matrix

    % EKF Loop
    for k = 1:N-1
        % Prediction Step
        x_pred = f(x_plus(:, k), u(:, k), dt); % Predict state
        F_k = F_func(x_plus(:, k), u(:, k), dt); % State transition Jacobian
        F_matrices(:, :, k) = F_k; % Store F_k for this step
        P_pred = F_k * P_plus(:, :, k) * F_k' + Q; % Predict covariance

        % Measurement Update Step
        H_k = H_func(x_pred); % Measurement Jacobian
        y_pred = h(x_pred); % Predicted measurement
        S = H_k * P_pred * H_k' + R; % Innovation covariance
        K = P_pred * H_k' / S; % Kalman gain

        % Innovation (residual)
        innovations(:, k+1) = y_meas(:, k+1) - y_pred;

        % Update state and covariance
        x_plus(:, k+1) = x_pred + K * innovations(:, k+1); % Updated state estimate
        P_plus(:, :, k+1) = (eye(state_dim) - K * H_k) * P_pred; % Updated covariance
    end

    % Store the final F matrix
    F_matrices(:, :, N) = F_func(x_plus(:, N), u(:, N), dt);
end
