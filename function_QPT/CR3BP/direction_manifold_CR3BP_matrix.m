function del_w_us_inter_full = direction_manifold_CR3BP_matrix(fin_qpos,p)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% calculate unstable and stable direction emanating from QPO

%%% input
%fin_qpos :invariant curve solution
%p :parameter dictionary

%%% output
%del_w_us_inter_full :unstable direction of manifolds
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% parameter dictionary
d = p("d");
N = p("N");
M = p("M");
mu = p("mu");

% options for ODE
options_ODE = odeset('RelTol',1e-13, 'AbsTol',1e-13);

% converged torus solution
T = fin_qpos(d*N*M+1);
rho = fin_qpos(d*N*M+2);

% fourier matrices
[FR, ~, ~] = fun_Fourier_matrix(rho,p);


% get unstable directions of full torus
% initialization
del_w_us_inter_full = zeros(d,N,M);
for i=1:M
    % invariant curve idx
    k=d*N*(i-1)+1:d*N*i;
    %    k=d*M*N*(i-1)+1:d*M*N*i;
    % invariant circle data
    X0 = fin_qpos(k);
    Mo = zeros(d*N,d*N);
    for j=1:N
        Y0 = [X0(d*j-(d-1):d*j);reshape(eye(d),[],1)];
        % ODE
        [~,Y] = ode113(@(t,x) fun_stm_cr3bp(t, x, mu),[0 T],Y0,options_ODE);
        % Monodromy Matrix of all states
        Mo(d*j-(d-1):d*j,d*j-(d-1):d*j) = reshape(Y(end,d+1:d+d^2),d,d);
    end
    % linearlization of stroboscopicmap
    Gw = FR*Mo;
    VD = zeros(d,N);
    for j=1:N
        [V,D] = eig(Gw(d*j-(d-1):d*j,d*j-(d-1):d*j));
        % eigenvalues
        D = diag(D);
        [~,idxmax] = max(abs(D));
        % eigenfunctions
        Vd = V(:,idxmax);
        if Vd(1)<0
            Vd = -Vd;
        end
        VD(:,j) = Vd;
    end
    del_w_us_inter_full(:,:,i) = VD;
end

end