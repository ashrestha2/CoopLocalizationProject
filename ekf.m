function [x_plus, P_plus, sigma, innovations, S, y_calc, F_matrices] = ekf(y_meas, x0, P0, Q, R, u, dt, N, f, h, F_func, H_func)
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

    I = eye(state_dim,state_dim);
    % EKF Loop
    for k = 1:N
        % Prediction Step
        % x_pred = f(x_plus(:, k), u(:, k), dt); % Predict state
        time_dist = [k-1 k]*dt;
        options = odeset('RelTol',1E-12,'AbsTol',1E-12);
        [t,X] = ode45(@(t,x) ugvEOM(t,x,u(:,k),zeros(6,1),0.5),time_dist,x_plus(:,k),options); 
        x_pred = X(end,:)';

        F_k = F_func(x_plus(:, k), u(:, k), dt); % State transition Jacobian
        F_matrices(:, :, k) = F_k; % Store F_k for this step
        P_pred = F_k * P_plus(:, :, k) * F_k' + Q; % Predict covariance

        % Measurement Update Step
        % H_k = H_func(x_pred); % Measurement Jacobian
        [~,~,H_k,~,~] = CT_to_DT(x_pred,0.5,2,12,-pi/18,pi/25,dt);

        y_pred = h(x_pred); % Predicted measurement
        S(:,:,k) = H_k * P_pred * H_k' + R; % Innovation covariance
        % check SK pos def
        try
        check = chol(H_k);
            disp('S is positive definite.');
        catch
            disp('S is NOT positive definite.');
        end

        % K = P_pred * H_k' / S(:,:,k); % Kalman gain
        % K = P_pred * H_k' * inv(S(:,:,k)); % Kalman gain
        K = (P_pred * H_k' * inv(S(:,:,k))); % Kalman gain

        % Innovation (residual)
        innovations(:, k) = y_meas(:, k) - y_pred;
        innovations(1,k) = wrapToPi(innovations(1,k));
        innovations(3,k) = wrapToPi(innovations(3,k));

        % Update state and covariance
        x_plus(:, k+1) = x_pred + K * innovations(:, k); % Updated state estimate
        P_plus(:, :, k+1) = (I - K * H_k) * P_pred; % Updated covariance
        sigma(:,k+1) = sqrt(diag(P_plus(:,:,k+1)));
        y_calc(:,k) = h(x_plus(:, k+1));

        % % check SK pos def
        % try
        % check = chol(P_plus(:, :, k+1));
        %     disp('S is positive definite.');
        % catch
        %     disp('S is NOT positive definite.');
        % end
    end

    % Store the final F matrix
    F_matrices(:, :, N) = F_func(x_plus(:, N), u(:, N), dt);
end
