function [F,G,H,M,Omega] = CT_to_DT(x_nom,L,vg,va,phi_g,omega_a,dt)

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

Gamma = eye(6,6);

% % Old code -- not sure if this does the same as the euler approx
% % Create continuous-time state-space model
% sys_ct = ss(A, B, C, D);
% 
% % Discretize the system
% sys_dt = c2d(sys_ct, dt);
% 
% % Extract discrete-time matrices
% F = sys_dt.A;
% G = sys_dt.B;
% H = sys_dt.C;
% M = sys_dt.D;

% use Euler approximation to calculate the DT matricies
F = eye(length(A)) + dt * A;
G = dt * B;
H = C;
M = D;
Omega = deltaT * Gamma;

% Display results
% disp('Discrete-time F matrix:');
% disp(F);
% disp('Discrete-time G matrix:');
% disp(G);
% disp('Discrete-time H matrix:');
% disp(H);
% disp('Discrete-time M matrix:');
% disp(M);

end