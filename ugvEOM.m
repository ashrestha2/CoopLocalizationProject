function [xdot] = ugvEOM(t,x,u,w,L)
    xi_g = x(1);
    eta_g = x(2);
    theta_g = x(3);
    xi_a = x(4);
    eta_a = x(5);
    theta_a = x(6);

    v_g = u(1);
    phi_g = u(2);
    v_a = u(3);
    w_a = u(4);

    w_xg = w(1);
    w_yg = w(2);
    w_wg = w(3);
    w_xa = w(4);
    w_ya = w(5);
    w_wa = w(6);

    xi_g_dot = v_g*cos(theta_g) + w_xg;
    eta_g_dot = v_g*sin(theta_g) + w_yg;
    theta_g_dot = (v_g/L)*tan(phi_g) + w_wg;

    xi_a_dot = v_a*cos(theta_a) + w_xa;
    eta_a_dot = v_a*sin(theta_a) + w_ya;
    theta_a_dot = w_a + w_wa;

    xdot = [xi_g_dot;eta_g_dot;theta_g_dot;xi_a_dot;eta_a_dot;theta_a_dot];
end