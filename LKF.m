function [x_plus, P_plus, innovation] = LKF(del_x0, P0, const, DT_mat_func, x_nom, y_nom, u_nom, y_meas, Q, R)

%%% Kalman Filter Function 
% Inputs: 
%   del_x0(n,1) = intial mean of the system 
%   P0(n,n) = initial covariance of the inital state
%   const = struct with all the constants 
%   DT_mat_func = function that gets DT matricies F_tilde, G_tilde,
%   H_tilde, M_tilde, Omega_tilde
%   R(p,p) = function call for a Measurement Error Covariance Matrix -- (p,p)
%   Q(n,n) = Process error covariance matrix (this is what should be tuned)
%   y_meas(p,T) = Measurement vector rows are all the meas and the columns 
%       are the samples
%   x_nom(n,T) = nomial states at each time 
%   y_nom(p,T) = nominal measurements at each time 

%
% Ouputs: 
%   x_plus(n,T) = a state matrix of all estimated states after they've been
%       filtered
%   P_plus(n,n,T) = a covariance 3d matrix with a stack of all the
%       covariance matricies after they've been filtered
%   innovation??/?
%   

    % initialize P(+) and del_x(+) 
    del_x_plus(:,1) = del_x0;
    P_plus(:,:,1) = P0;
    
    % determine number of states and samples
    n = length(del_x_plus); % how many states
    T = length(y_meas); % how many samples were taken

    %create some matricies to be used
    I = eye(n);

    % running through the loop for every time step 
    for k = 1:T %k = time step
        [F_tilde,G_tilde,H_tilde,M_tilde,omega_tilde] = DT_matricies_func(x_nom,const.L,const.vg0,const.va0,const.phi_g0,const.omega_a0,const.deltaT);

        %%%%%%%%%%%%%%%%%%%%%%%%
        %%% prediction step section 
        del_u(:,k) = u(:,k+1) - u_nom(:,k+1); %WHERE THE HECK DOES UK+1 COME FROM -- 0???
        del_x_minus(:,k+1) = F_tilde * del_x_plus(:,k) + G_tilde * del_u(:,k);
        P_minus = F_tilde * P_plus(:,:,k) * F_tilde' + omega_tilde * Q * omega_tilde';
        innovation(:,:,k) = H_tilde*P_minus*H_tilde'+R;
        K(:,:,k+1) = P_minus * H_tilde' * inv(innovation(:,:,k));
    
        % % transition from minus to plus (testing the prediciton step)
        % del_x_plus(:,k+1) = del_x_minus(:,k+1);
        % P_plus(:,:,k+1) = P_minus(:,:,k+1);
    
        %%%%%%%%%%%%%%%%%%%%%%%
        %%% correction step 
        del_y_meas(:,k) = y_meas(:,k) - y_nom(:,k); % y_meas given to us in mat file y_nom thru findYnom func with x_nom
        del_x_plus(:,k+1) = del_x_minus(:,k+1) + K(:,:,k+1) * (del_y_meas(:,k) - (H_tilde * del_x_minus(:,k+1)));
        P_plus(:,:,k+1) = (I - K(:,:,k+1) * H_tilde) * P_minus;
    
    end

end