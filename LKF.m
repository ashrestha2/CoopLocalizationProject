function [x_plus, P_plus, innovation] = Kalman_Linearized_Filter(x0, P0, ...
    F_tilde, G_tilde, del_u, Q, omega_tilde, y_meas, x_nom, u_nom, y_nom)
% how to get all the nom -- should push in functions????


%%% Kalman Filter Function 
% Inputs: 
%   u0(n,1) = intial mean of the system 
%   P0(n,n,1) = initial covariance of the inital state

% SHOULD I MAKE THEM ALL FUNCTION CALLS?
%   F_tilde = function call for a DT transform matrix from one 
%       perturbation state to the next perturbation state -- (n,n) matrix
%   G_tilde = function call for a DT transform matrix from one input to 
%       the next perturbation state -- (n,m) matrix
%   H = function call or a DT transform matrix from states to
%       measurements -- (p,n) matrix
%   R = function call for a Measurement Error Covariance Matrix -- (p,p)
%   Q(n,n) = Process error covariance matrix (this is what should be tuned)
%   Omega_tilde() = ____________
%   y_meas(p,T) = Measurement vector rows are all the meas and the columns 
%       are the samples
%   x_nom(n,T) = nomial states at each time 
%
% Ouputs: 
%   x_plus(n,T) = a state matrix of all estimated states after they've been
%       filtered
%   P_plus(n,n,T) = a covariance 3d matrix with a stack of all the
%       covariance matricies after they've been filtered
%   innovation??/?
%   

% initialize P(+) and del_x(+) 
    del_x_plus(:,1) = x0;
    P_plus(:,:,1) = P0;
    
    % determine number of states and samples
    n = length(del_x_plus); % how many states
    T = length(y_meas); % how many samples were taken

    %create some matricies to be used
    I = eye(n);


    % running through the loop for every time step 
    for k = 1:T %k = time step
    
        %%%%%%%%%%%%%%%%%%%%%%%%
        %%% prediction step section 
        del_x_minus(:,k+1) = F_tilde * del_x_plus(:,k) + G_tilde * del_u(:,k);
        P_minus(:,:,k+1) = F_tilde * P_plus(:,:,k) * F_tilde' + omega_tilde * Q * omega_tilde';
        innovation(:,:,k) = H_tilde*P_minus(:,:,k+1)*H_tilde'+R;
        K(:,:,k+1) = P_minus(:,:,k+1) * H_tilde' * inv(innovation(:,:,k));
    
        del_u(:,k) = u(:,k+1) - u_nom(:,k+1); %WHERE THE HECK DOES UK+1 COME FROM
        % % transition from minus to plus (testing the prediciton step)
        % del_x_plus(:,k+1) = del_x_minus(:,k+1);
        % P_plus(:,:,k+1) = P_minus(:,:,k+1);
    
        %%%%%%%%%%%%%%%%%%%%%%%
        %%% correction step 
        del_x_plus(:,k+1) = del_x_minus(:,k+1) + K(:,:,k+1) * (y_meas(:,k) - H_tilde * del_x_minus(:,k+1));
        P_plus(:,:,k+1) = (I - K(:,:,k+1) * H_tilde) * P_minus(:,:,k+1);
        del_y(:,k+1) = y_meas(:,k+1) - y_nom; 
    
    end

endfunction [x_plus, P_plus, innovation] = Kalman_Linearized_Filter(x0, P0, ...
    F_tilde, G_tilde, del_u, Q, omega_tilde, y_meas, x_nom, u_nom, y_nom)
% how to get all the nom -- should push in functions????


%%% Kalman Filter Function 
% Inputs: 
%   u0(n,1) = intial mean of the system 
%   P0(n,n,1) = initial covariance of the inital state

% SHOULD I MAKE THEM ALL FUNCTION CALLS?
%   F_tilde = function call for a DT transform matrix from one 
%       perturbation state to the next perturbation state -- (n,n) matrix
%   G_tilde = function call for a DT transform matrix from one input to 
%       the next perturbation state -- (n,m) matrix
%   H = function call or a DT transform matrix from states to
%       measurements -- (p,n) matrix
%   R = function call for a Measurement Error Covariance Matrix -- (p,p)
%   Q(n,n) = Process error covariance matrix (this is what should be tuned)
%   Omega_tilde() = ____________
%   y_meas(p,T) = Measurement vector rows are all the meas and the columns 
%       are the samples
%   x_nom(n,T) = nomial states at each time 
%
% Ouputs: 
%   x_plus(n,T) = a state matrix of all estimated states after they've been
%       filtered
%   P_plus(n,n,T) = a covariance 3d matrix with a stack of all the
%       covariance matricies after they've been filtered
%   innovation??/?
%   

% initialize P(+) and del_x(+) 
    del_x_plus(:,1) = x0;
    P_plus(:,:,1) = P0;
    
    % determine number of states and samples
    n = length(del_x_plus); % how many states
    T = length(y_meas); % how many samples were taken

    %create some matricies to be used
    I = eye(n);


    % running through the loop for every time step 
    for k = 1:T %k = time step
    
        %%%%%%%%%%%%%%%%%%%%%%%%
        %%% prediction step section 
        del_x_minus(:,k+1) = F_tilde * del_x_plus(:,k) + G_tilde * del_u(:,k);
        P_minus(:,:,k+1) = F_tilde * P_plus(:,:,k) * F_tilde' + omega_tilde * Q * omega_tilde';
        innovation(:,:,k) = H_tilde*P_minus(:,:,k+1)*H_tilde'+R;
        K(:,:,k+1) = P_minus(:,:,k+1) * H_tilde' * inv(innovation(:,:,k));
    
        del_u(:,k) = u(:,k+1) - u_nom(:,k+1); %WHERE THE HECK DOES UK+1 COME FROM
        % % transition from minus to plus (testing the prediciton step)
        % del_x_plus(:,k+1) = del_x_minus(:,k+1);
        % P_plus(:,:,k+1) = P_minus(:,:,k+1);
    
        %%%%%%%%%%%%%%%%%%%%%%%
        %%% correction step 
        del_x_plus(:,k+1) = del_x_minus(:,k+1) + K(:,:,k+1) * (y_meas(:,k) - H_tilde * del_x_minus(:,k+1));
        P_plus(:,:,k+1) = (I - K(:,:,k+1) * H_tilde) * P_minus(:,:,k+1);
        del_y(:,k+1) = y_meas(:,k+1) - y_nom; 
    
    end

endfunction [x_plus, P_plus, innovation] = Kalman_Linearized_Filter(x0, P0, ...
    F_tilde, G_tilde, del_u, Q, omega_tilde, y_meas, x_nom, u_nom, y_nom)
% how to get all the nom -- should push in functions????


%%% Kalman Filter Function 
% Inputs: 
%   u0(n,1) = intial mean of the system 
%   P0(n,n,1) = initial covariance of the inital state

% SHOULD I MAKE THEM ALL FUNCTION CALLS?
%   F_tilde = function call for a DT transform matrix from one 
%       perturbation state to the next perturbation state -- (n,n) matrix
%   G_tilde = function call for a DT transform matrix from one input to 
%       the next perturbation state -- (n,m) matrix
%   H = function call or a DT transform matrix from states to
%       measurements -- (p,n) matrix
%   R = function call for a Measurement Error Covariance Matrix -- (p,p)
%   Q(n,n) = Process error covariance matrix (this is what should be tuned)
%   Omega_tilde() = ____________
%   y_meas(p,T) = Measurement vector rows are all the meas and the columns 
%       are the samples
%   x_nom(n,T) = nomial states at each time 
%
% Ouputs: 
%   x_plus(n,T) = a state matrix of all estimated states after they've been
%       filtered
%   P_plus(n,n,T) = a covariance 3d matrix with a stack of all the
%       covariance matricies after they've been filtered
%   innovation??/?
%   

% initialize P(+) and del_x(+) 
    del_x_plus(:,1) = x0;
    P_plus(:,:,1) = P0;
    
    % determine number of states and samples
    n = length(del_x_plus); % how many states
    T = length(y_meas); % how many samples were taken

    %create some matricies to be used
    I = eye(n);


    % running through the loop for every time step 
    for k = 1:T %k = time step
    
        %%%%%%%%%%%%%%%%%%%%%%%%
        %%% prediction step section 
        del_x_minus(:,k+1) = F_tilde * del_x_plus(:,k) + G_tilde * del_u(:,k);
        P_minus(:,:,k+1) = F_tilde * P_plus(:,:,k) * F_tilde' + omega_tilde * Q * omega_tilde';
        innovation(:,:,k) = H_tilde*P_minus(:,:,k+1)*H_tilde'+R;
        K(:,:,k+1) = P_minus(:,:,k+1) * H_tilde' * inv(innovation(:,:,k));
    
        del_u(:,k) = u(:,k+1) - u_nom(:,k+1); %WHERE THE HECK DOES UK+1 COME FROM
        % % transition from minus to plus (testing the prediciton step)
        % del_x_plus(:,k+1) = del_x_minus(:,k+1);
        % P_plus(:,:,k+1) = P_minus(:,:,k+1);
    
        %%%%%%%%%%%%%%%%%%%%%%%
        %%% correction step 
        del_x_plus(:,k+1) = del_x_minus(:,k+1) + K(:,:,k+1) * (y_meas(:,k) - H_tilde * del_x_minus(:,k+1));
        P_plus(:,:,k+1) = (I - K(:,:,k+1) * H_tilde) * P_minus(:,:,k+1);
        del_y(:,k+1) = y_meas(:,k+1) - y_nom; 
    
    end

end