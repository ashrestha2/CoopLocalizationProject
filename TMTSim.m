function [time_iter,x_array,y] = TMTSim(const,Qtrue,Rtrue)
    t_int = 0:const.deltaT:100;
    u = [const.v_g0;const.phi_g0;const.v_a0;const.w_a0];
    perturb_x0 = [0;1;0;0;0;0.1];
    x_array(1,:) = const.x0 + perturb_x0;
    time_iter(1) = 0;
    for i = 1: length(t_int) - 1
        dt = [t_int(i) t_int(i+1)];
        n = length(const.x0);
        w = mvnrnd(zeros(n,1),Qtrue)';
        options = odeset('RelTol',1E-12,'AbsTol',1E-12);
        [t,X] = ode45(@(t,x) ugvEOM(t,x,u,w,const.L),dt,x_array(i,:),options); 
        time_iter(i+1) = t(end);
        x_array(i+1,:) = X(end,:);
    end
    for i = 2:length(t_int)
        p = 5;
        v = mvnrnd(zeros(p,1),Rtrue)';
        y(i-1,:) = findYnom(x_array(i,:)) + v;
    end
    y = y';
end