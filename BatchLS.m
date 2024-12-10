function [x_bls, P_bls] = BatchLS(y_meas_start,R,x_nom,const,DT_mat_func)

% C batch weighted 
numstates = 5;
time_steps = 10;

%unwrap y
T = length(y_meas_start);
y_new = reshape(y_meas_start,[T*numstates 1]);

bigH = [];
bigR = [];

for j = 1:time_steps
    [~,~,H,~,~] = DT_mat_func(x_nom(j,:),const.L,const.v_g0,const.v_a0,const.phi_g0,const.w_a0,const.deltaT);
    bigH = [bigH;eye(numstates)];
    bigR = blkdiag(bigR,R);
end

% find batch least squares
x_bls(:,1) = inv(bigH'*inv(bigR)*bigH)*(bigH'*inv(bigR))*y_new(1:numstates*time_steps,:)

% error covariance matrix
P_bls(:,:)=inv(bigH'*inv(bigR)*bigH);

end