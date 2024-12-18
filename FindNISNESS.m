function [epsNEESbar,r1x,r2x,epsNISbar,r1y,r2y, NEES, NIS] = FindNISNESS(N,del_x0,P0,x_nom,y_nom,DT_mat_func,const,Qtrue,Rtrue,endTime)
    n = 6;
    p = 5;
    for i = 1:N
        % Generate true trajectory and measurements from system
        [t,xtrue,ytrue] = TMTSim(const,Qtrue,Rtrue,endTime);

        % Kalman filter 
        [x_plus, P_plus, Sk, y_minus, sigma, innovation] = LKF(del_x0, P0, const, DT_mat_func, x_nom, y_nom, ytrue, Qtrue, Rtrue);


        % Calculate NEES and NIS
        for j = 1: length(t)
            error_x(:,i,j) = xtrue(j,:) - x_plus(:,j)';
            error_x(3,i,j) = wrapToPi(error_x(3,i,j));
            error_x(6,i,j) = wrapToPi(error_x(6,i,j));
            NEES(i,j) = (error_x(:,i,j)')*inv(P_plus(:,:,j))*(error_x(:,i,j));
            if j > 2
                NIS(i,j-1) = (innovation(:,j-1))'*inv(Sk(:,:,j-1))*(innovation(:,j-1));
            end
        end
    end

    % NEES Test:
    epsNEESbar = mean(NEES,1);
    alphaNEES = 0.05;
    Nnx = N*n;
    r1x = chi2inv(alphaNEES/2, Nnx)./N;
    r2x = chi2inv(1-alphaNEES/2, Nnx)./N;
   
    figure(15);
    plot(epsNEESbar,'ro','MarkerSize',6,'LineWidth',2),hold on
    plot(r1x*ones(size(epsNEESbar)),'r--','LineWidth',2)
    plot(r2x*ones(size(epsNEESbar)),'r--','LineWidth',2)
    ylabel('NEES statistic, $\bar{\epsilon}_x$','FontSize',14,'Interpreter','latex')
    xlabel('time step, k','FontSize',14)
    title('NEES Estimation Results','FontSize',14)
    legend('NEES @ time k', 'r_1 bound', 'r_2 bound')
    ylim([0 15])
    
    epsNISbar = mean(NIS,1);
    alphaNIS = 0.1;
    Nny = N*p;
    r1y = chi2inv(alphaNIS/2,Nny)./N;
    r2y = chi2inv(1-alphaNIS/2,Nny)./N;

    figure(16);
    plot(epsNISbar,'bo','MarkerSize',6,'LineWidth',2),hold on
    plot(r1y*ones(size(epsNISbar)),'b--','LineWidth',2)
    plot(r2y*ones(size(epsNISbar)),'b--','LineWidth',2)
    ylabel('NIS statistic, $\bar{\epsilon}_y$','FontSize',14,'Interpreter','latex')
    xlabel('time step, k','FontSize',14)
    title('NIS Estimation Results','FontSize',14)
    legend('NIS @ time k', 'r_1 bound', 'r_2 bound')
    ylim([0 15])


figure(17);
subplot(5,1,1); hold on;
plot(innovation(1,:));
plot(+2*sqrt(reshape(Sk(1,1,2:end),[1,length(Sk)-1])),'b--',LineWidth=1.2)
plot(-2*sqrt(reshape(Sk(1,1,2:end),[1,length(Sk)-1])),'b--',LineWidth=1.2)
subplot(5,1,2); hold on;
plot(innovation(2,:));
plot(+2*sqrt(reshape(Sk(2,2,2:end),[1,length(Sk)-1])),'b--',LineWidth=1.2)
plot(-2*sqrt(reshape(Sk(2,2,2:end),[1,length(Sk)-1])),'b--',LineWidth=1.2)
subplot(5,1,3); hold on;
plot(innovation(3,:));
plot(+2*sqrt(reshape(Sk(3,3,2:end),[1,length(Sk)-1])),'b--',LineWidth=1.2)
plot(-2*sqrt(reshape(Sk(3,3,2:end),[1,length(Sk)-1])),'b--',LineWidth=1.2)
subplot(5,1,4); hold on;
plot(innovation(4,:));
plot(+2*sqrt(reshape(Sk(4,4,2:end),[1,length(Sk)-1])),'b--',LineWidth=1.2)
plot(-2*sqrt(reshape(Sk(4,4,2:end),[1,length(Sk)-1])),'b--',LineWidth=1.2)
subplot(5,1,5); hold on;
plot(innovation(5,:));
sgtitle('innovation')


end