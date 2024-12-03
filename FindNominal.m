function x_dot = FindNominal(t,x,v_g0, v_a0, L, phi_g0, w_a0)
% pull out states
theta_g = x(3);
theta_a = x(6);


% xdot = f(x,const)
nom_east_g = v_g0*cos(theta_g);
nom_north_g = v_g0*sin(theta_g);
nom_theta_g = wrapToPi(v_g0/L * tan(phi_g0)) ;
nom_east_a = v_a0*cos(theta_a);
nom_north_a = v_a0*sin(theta_a);
nom_theta_a = wrapToPi(w_a0) ;

x_dot = [nom_east_g; nom_north_g; nom_theta_g; nom_east_a; nom_north_a; nom_theta_a];

end