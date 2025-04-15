function del_w_us_inter_full = direction_manifold(fin_qpos,p,FR)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% calculate unstable and stable direction emanating from QPO
%% By:Soi Yamaguchi

%%% input
%fin_qpos :invariant curve solution
%p :parameter dictionary
%FR :Fourier Matirx dictionary

%%% output
%del_w_us_inter_full :unstable direction of manifolds
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% parameter dictionary
d = p("d");
N = p("N");
M = p("M");
mu = p("mu");

% options for ODE
options = odeset('RelTol',1e-13, 'AbsTol',1e-13);

% converged torus solution
% T
T = fin_qpos(d*N*M+1);
% rho
rho = fin_qpos(d*N*M+2);

% fourier matrices
Fr = FR{1,2};
IFr = FR{2,2};
Q = FR{4,2};
Q = rotation_matrix(Q,rho,p);
R = IFr*Q*Fr;

% get unstable directions of full torus
% initialization
del_w_us_inter_full = zeros(d,N,M);
for i=1:M
    % invariant curve idx
    k=d*N*(i-1)+1:d*N*i;
    % invariant circle data
    X0 = fin_qpos(k);
    Mo = zeros(d*N,d*N);
    for j=1:N
        Y0 = [X0(d*j-5:d*j);reshape(eye(d),[],1)];
        % ODE
        [~,Y] = ode113(@(t,x) fun_stm_cr3bp(t,x,mu),[0 T],Y0,options);
        % Monodromy Matrix of all states
        Mo(d*j-(d-1):d*j,d*j-(d-1):d*j) = reshape(Y(end,d+1:d+d^2),d,d);
    end
    % linearlization of stroboscopicmap
    Gw = R*Mo;
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