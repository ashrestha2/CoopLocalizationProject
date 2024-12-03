function y_nom = findYnom(x_nom)
    % pull out all the states
    e_g(:,1) = x_nom(:,1);
    n_g(:,1) = x_nom(:,2);
    theta_g(:,1) = x_nom(:,3);
    e_a(:,1) = x_nom(:,4);
    n_a(:,1) = x_nom(:,5);
    theta_a(:,1) = x_nom(:,6);
    
    % calculate each row of the measurements
    y_nom(1,:) = atan2((n_a-n_g),(e_a-e_g)) - theta_g;
    y_nom(2,:) = sqrt((e_g-e_a).^2+(n_g-n_a).^2);
    y_nom(3,:) = atan2((n_g-n_a),(e_g-e_a)) - theta_a;
    y_nom(4,:) = e_a;
    y_nom(5,:) = n_a;

end