function del_w_us_interpolation = direction_manifold_interpolation(fin_qpos, p)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate the direction manifold interpolation
% By: Soichiro Shin

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%% input
% fin_qpos :invariant curve solution
% p :parameter

% %%% output
% del_w_us_interpolation :unstable direction manifold interpolation

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% parameter dictionary
d = p("d");
N = p("N");
M = p("M");
mu = p("mu");
num_iter = p("num_iter");


% option for ODE
options_ODE = odeset('RelTol',1e-13,'AbsTol',1e-13);

% converged torus solution
T = fin_qpos(d*N*M+1);
rho = fin_qpos(d*N*M+2);

% fourier matrices
FR = fun_Fourier_interpolation(rho,p);

% get unstabel directions of full torus
% initialization
del_w_us_interpolation = zeros(d,num_iter,num_iter);
for i = 1:num_iter
    % invariant curve idx
    k = d*num_iter*(i-1)+1:d*num_iter*i;
    % invariant circle data
    X0 = fin_qpos(k);
    Mo = zeros(d*num_iter,d*num_iter);
    for j = 1:num_iter
        Y0 = [X0(d*j-(d-1):d*j); reshape(eye(d),[],1)];
        % ODE
        [~,Y] = ode113(@(t,x) fun_stm_cr3bp(t,x,mu),[0 T],Y0,options_ODE);
        % Monodromy matrix of all states
        Mo(d*j-(d-1):d*j,d*j-(d-1):d*j) = reshape(Y(end,d+1:d+d^2),d,d);
    end
    % lineralization of stroboscopic map
    Gw = FR*Mo;
    VD = zeros(d,num_iter);
    for j = 1:num_iter
        % eigenvalue and eigenvector
        [V,D] = eig(Gw(d*j-(d-1):d*j,d*j-(d-1):d*j));
        % eigenvalue
        D = diag(D);
        [~,idmax] = max(abs(D));
        % eigenfunctions
        Vd = V(:,idmax);
        if Vd(1) < 0
            Vd = -Vd;
        end
        VD(:,j) = Vd;
    end
    del_w_us_interpolation(:,:,i) = VD;
end

end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
