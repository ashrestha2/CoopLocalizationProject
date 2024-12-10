function [x_plus, P_plus, innovation] = Kalman_Filter(u0, P0, F, H, R, Q, y_meas)
%%% Kalman Filter Function 
% Inputs: 
%   u0(n,1) = intial mean of the system 
%   P0(n,n,1) = initial covariance of the inital state
%   F(n,n) = DT transform matrix from one state to the next state 
%       F(n,n,T) = stack of DT transform matrcies from one state to the
%           next specific to each time k
%   H(p,n) = transform matrix from states to measurements
%       H(p,n,T) stacked transformation matricies from states to each
%           different k measurement 
%   R(p,p) = Measurement Error Covariance Matrix 
%       R(p,p,T) = block matrix of meas error covariance matricies for each
%           different k measurement
%   Q(n,n) = Process error covariance matrix (this is what should be tuned)
%   y_meas(p,T) = Measurement vector rows are all the meas and the columns 
%       are the samples
%
% Ouputs: 
%   x_plus(n,T) = a state matrix of all estimated states after they've been
%       filtered
%   P_plus(n,n,T) = a covariance 3d matrix with a stack of all the
%       covariance matricies after they've been filtered
%   innovation??/?
%   

    % initialize
    x_plus(:,1) = u0;
    P_plus(:,:,1) = P0;
    
    % determine number of states and samples
    n = length(u0); % how many states
    T = length(y_meas); % how many measuremetns per sample

    %create some matricies to be used
    I = eye(n);

    %%% Make sure all the F, H and R matricies are all 3d matricies even if
    %%% they don't change with time
    % check of F changes with time -- if not create a 3d matrix so the for
    % loop can handle it 
    if  length(size(F)) == 2
        F = repmat(F,[1,1,T]);
    end
    % check of H changes with time -- if not create a 3d matrix so the for
    % loop can handle it 
    if  length(size(H)) == 2
        H = repmat(H,[1,1,T]);
    end
    % check of R changes with time -- if not create a 3d matrix so the for
    % loop can handle it 
    if length(size(R)) == 2
        R = repmat(R,[1,1,T]);
    end
    
%%% NEED TO CHECK IF INDEXING H R CORRECTLY SINCE IF I INDEX H OR R
%%% STARTING AT K+1 IT DOESN'T WORK 

    % running through the loop
    for k = 1:T
    
        %%%%%%%%%%%%%%%%%%%%%%%%
        %%% prediction step section 
        x_minus(:,k+1) = F(:,:,k) * x_plus(:,k);
        P_minus(:,:,k+1) = F(:,:,k) * P_plus(:,:,k) * F(:,:,k)' + Q;
        innovation(:,:,k) = H(:,:,k)*P_minus(:,:,k+1)*H(:,:,k)'+R(:,:,k);
        K(:,:,k+1) = P_minus(:,:,k+1) * H(:,:,k)' * inv(innovation(:,:,k));
    
        % % transition from minus to plus (testing the prediciton step)
        % x_plus(:,k+1) = x_minus(:,k+1);
        % P_plus(:,:,k+1) = P_minus(:,:,k+1);
    
        %%%%%%%%%%%%%%%%%%%%%%%
        %%% correction step 
        x_plus(:,k+1) = x_minus(:,k+1) + K(:,:,k+1) * (y_meas(:,k) - H(:,:,k) * x_minus(:,k+1));
        P_plus(:,:,k+1) = (I - K(:,:,k+1) * H(:,:,k)) * P_minus(:,:,k+1);
    
    end
end