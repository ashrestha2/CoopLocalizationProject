% NOminal System Parameters
L = 0.5; % Vehicle length (m)
vg = 2; % Ground vehicle speed (m/s)
va = 12; % Aerial vehicle speed (m/s)
phi_g = -pi/18; % Steering angle (rad)
omega_a = pi/25; % Angular rate (rad/s)
dt = 0.1; % Sampling time (s)

% Nominal trajectory
x_nom = [10; 0; pi/2; -60; 0; -pi/2]; % [zeta_g (1), eta_g (2), theta_g (3), zeta_a (4), eta_a (5), theta_a (6)]

% ABCD CT matrices 
A = [
    0, 0, -vg*sin(x_nom(3)), 0, 0, 0;
    0, 0, vg*cos(x_nom(3)), 0, 0, 0;
    0, 0, 0, 0, 0, 0;
    0, 0, 0, 0, 0, -va*sin(x_nom(6));
    0, 0, 0, 0, 0, va*cos(x_nom(6));
    0, 0, 0, 0, 0, 0
];

B = [
    cos(x_nom(3)), 0, 0, 0;
    sin(x_nom(3)), 0, 0, 0;
    tan(phi_g)/L, vg*(1 + tan(phi_g)^2)/L, 0, 0;
    0, 0, cos(x_nom(6)), 0;
    0, 0, sin(x_nom(6)), 0;
    0, 0, 0, 1
];

C = [(x_nom(5) - x_nom(2))/((x_nom(4) - x_nom(1))^2 + (x_nom(5) - x_nom(2))^2),     -(x_nom(4) - x_nom(1))/((x_nom(4) - x_nom(1))^2 + (x_nom(5) - x_nom(2))^2),     -1, -(x_nom(5) - x_nom(2))/((x_nom(4) - x_nom(1))^2 + (x_nom(5) - x_nom(2))^2),     (x_nom(4) - x_nom(1))/((x_nom(4) - x_nom(1))^2 + (x_nom(5) - x_nom(2))^2),      0;
    (x_nom(1) - x_nom(4))/sqrt((x_nom(1) - x_nom(4))^2 + (x_nom(2) - x_nom(5))^2),  (x_nom(2) - x_nom(5))/sqrt((x_nom(1) - x_nom(4))^2 + (x_nom(2) - x_nom(5))^2),   0, -(x_nom(1) - x_nom(4))/sqrt((x_nom(1) - x_nom(4))^2 + (x_nom(2) - x_nom(5))^2), -(x_nom(2) - x_nom(5))/sqrt((x_nom(1) - x_nom(4))^2 + (x_nom(2) - x_nom(5))^2), 0;
    -(x_nom(2) - x_nom(5))/((x_nom(1) - x_nom(4))^2 + (x_nom(2) - x_nom(5))^2),     (x_nom(1) - x_nom(4))/((x_nom(1) - x_nom(4))^2 + (x_nom(2) - x_nom(5))^2),       0, (x_nom(2) - x_nom(5))/((x_nom(1) - x_nom(4))^2 + (x_nom(2) - x_nom(5))^2),      -(x_nom(1) - x_nom(4))/((x_nom(1) - x_nom(4))^2 + (x_nom(2) - x_nom(5))^2),     -1;
    0, 0, 0, 1, 0, 0;
    0, 0, 0, 0, 1, 0
];

D = zeros(size(C, 1), size(B, 2)); 

% Create continuous-time state-space model
sys_ct = ss(A, B, C, D);

% Discretize the system
sys_dt = c2d(sys_ct, dt);

% Extract discrete-time matrices
F = sys_dt.A;
G = sys_dt.B;
H = sys_dt.C;
M = sys_dt.D;

% Display results
disp('Discrete-time F matrix:');
disp(Ad);
disp('Discrete-time G matrix:');
disp(Bd);
disp('Discrete-time H matrix:');
disp(Cd);
disp('Discrete-time M matrix:');
disp(Dd);
