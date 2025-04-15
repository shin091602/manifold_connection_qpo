function [hinvcir, check_accuracy_pink] = plot_invariant_curve_CR3BP_matrix(x_star_pert,zpo,z0,p)
%% plot invariatnt curve
%%% input
%zpo :periodic orbit
%z0 :initial invariant curve solution
%p :parameter

%%% output
%hinvcir :initial invariant curve plot

%% DICTIONARY OPEN
d = p("d");
N = p("N");
M = p("M");
mu = p("mu");

%% OPTIONS ODE
options_ODE = odeset('RelTol',1e-13, 'AbsTol',1e-13);

% invariant curve
% X0 = z0(1:d*N*M);
% X0 = reshape(X0, d, N, M); % =m
X0 = x_star_pert; %(d,N,M)
T = z0(d*N*M+1);
rho = z0(d*N*M+2);

%% fourier matrices
[FR, ~, ~] = fun_Fourier_matrix(rho, p);

% time span
inti = linspace(0,T,M+1);
inti = [inti,T]; %いらなくね

%% make orbit
Y_end = zeros(d,N,M);
X1 = zeros(d*N, M);

for i = 1:M
    for j = 1:N
        % time span
        ts = [inti(i) inti(i+1)];
        % propagete orbit for 1T = blue
        [~,Y] = ode113(@(t,x) fun_cr3bp(t,x,mu),ts,X0(:,j,i),options_ODE);
        Y = Y';
        Y_end(:,j,i) = Y(:,end);
    end
end

Yf = reshape(Y_end, d*N, M);

for i = 1:M
    % DFT → rotate → IDFT = make green
    X1(:,i) = FR*Yf(:,i);
end

%% plot (i=Mのもののみ)
if p("HLflg") %Halo
    hinvcir = figure();
    hold on
    grid on
    axis equal
    view(10, 10)
    xlabel('$x$[-]');
    ylabel('$y$[-]');
    zlabel('$z$[-]');
    % periodic orbit
    scatter3(zpo(1),zpo(2),zpo(3),400,'r.');
    % invariant circle
    for j=1:p("N")
        % for j=1:3
        finvcir = scatter3(X0(1,j,1),X0(2,j,1),X0(3,j,1), 100,'m.');
        fstrobo = scatter3(Y_end(1,j,end),Y_end(2,j,end), Y_end(3,j,end), 100,'b.');
        fR = scatter3(X1(d*j-(d-1),end),X1(d*j-(d-2),end),X1(d*j-(d-3),end),100,'g.');
    end
    legend([finvcir,fstrobo,fR],{'$\Psi(\theta)$','$\phi_{T}(\Psi(\theta))$','$R_{-\rho}[\phi_{T}(\Psi(\theta))]$'},'Location','northoutside');
    lgd = legend;
    lgd.NumColumns = 3;
    hold off

else %Lyapunov
    if p("rho_type") %in-plane
        hinvcir = figure();
        hold on
        grid on
        axis equal
        xlabel('$x$[-]');
        ylabel('$y$[-]');
        % periodic orbit
        scatter(zpo(1),zpo(2),400,'r.');

        % invariant circle
        for j=1:p("N")
            finvcir = scatter(X0(1,j,1),X0(2,j,1), 100,'m.');
            fstrobo = scatter(Y_end(1,j,end),Y_end(2,j,end), 100,'b.');
            fR = scatter(X1(d*j-(d-1),end),X1(d*j-(d-2),end),100,'g.');
            % finvcir = scatter(X(d*j-(d-3)),X(d*j),100,'m.');
            % fstrobo = scatter(Yf(d*j-(d-3),end),Yf(d*j,end),100,'b.');
            % fR = scatter(RYf(d*j-(d-3),end),RYf(d*j,end),100,'g.');
        end
        legend([finvcir,fstrobo,fR],{'$\Psi(\theta)$','$\phi_{T}(\Psi(\theta))$','$R_{-\rho}[\phi_{T}(\Psi(\theta))]$'},'Location','northoutside');
        lgd = legend;
        lgd.NumColumns = 3;
        hold off
    else % out-plane
        hinvcir = figure();
        hold on
        grid on
        axis equal
        xlabel('$z$[-]');
        ylabel('$\dot{z}$[-]');

        % periodic orbit
        scatter(zpo(3),zpo(6),400,'r.');

        % invariant circle
        for j=1:p("N")
            finvcir = scatter(X0(3,j,1),X0(6,j,1), 100,'m.');
            fstrobo = scatter(Y_end(3,j,end),Y_end(6,j,end), 100,'b.');
            fR = scatter(X1(d*j-(d-3),end),X1(d*j-(d-6),end),100,'g.');
            % finvcir = scatter(X(d*j-(d-3)),X(d*j),100,'m.');
            % fstrobo = scatter(Yf(d*j-(d-3),end),Yf(d*j,end),100,'b.');
            % fR = scatter(RYf(d*j-(d-3),end),RYf(d*j,end),100,'g.');
        end
        legend([finvcir,fstrobo,fR],{'$\Psi(\theta)$','$\phi_{T}(\Psi(\theta))$','$R_{-\rho}[\phi_{T}(\Psi(\theta))]$'},'Location','northoutside');
        lgd = legend;
        lgd.NumColumns = 3;
        hold off
    end
end

X0_vec = reshape(X0, d*N,M);
X1_vec = reshape(X1, d*N,M);
check_accuracy_pink = norm(X0_vec(:,1) - X1_vec(:,end));
end