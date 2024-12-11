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

% Define parameters
y_meas = randn(5, 100); % Measurements (5 x 100)
x0 = [10; 0; pi/2; -60; 0; -pi/2]; % Initial state
P0 = eye(6); % Initial covariance
Q = 0.01 * eye(6); % Process noise covariance
R = 0.1 * eye(5); % Measurement noise covariance
u = repmat([2; -pi/18; 12; pi/25], 1, 100); % Constant control inputs
dt = 0.1; % Sampling time
N = 100; % Number of time steps

% Call EKF
[x_plus, P_plus, innovations] = ekf(y_meas, x0, P0, Q, R, u, dt, N, f, h, F_func, H_func);
