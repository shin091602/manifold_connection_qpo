function [z0,phi0,Ud0,rho,Ec] = fun_center_manifold_cr3bp(zpo,p,FR)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% calculate invariant circle and rotation number
%% By:Soi Yamaguchi

%%% input
%zpo :initial data of periodidxcen orbit
%p :parameter didxcentionary
%FR :fourier matrix dictionary

%%% output
%z0 :initial guess of invariant circle
%phi0 :initial tangent vector
%Ud0 :torus function dictionary
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% DICTIONARY
d = p("d");
M = p("M");
N = p("N");
K = p("K");
mu = p("mu");

%% OPTIONS ODE
options = odeset('RelTol',1e-13, 'AbsTol',1e-13);

%% calculate eigenvalues of the monodromy matrix
% Propogate orbit for one period and gather data at N points
% periodidxcen orbit state
xpo = zpo(1:d);
% periodidxcen orbit period
Tpo = zpo(end);
% time vector
tk = linspace(0,Tpo,M+1);
tk = tk(1,1:end-1);
% invariant circle angles
thtj = linspace(0,2*pi,N+1);
thtj = thtj(1,1:end-1);

% initial state including initial STM
XP0 = [xpo; reshape(eye(d), [], 1)];
[~,Y] = ode113(@(t,x) fun_stm_cr3bp(t,x,mu),[0 Tpo],XP0,options);

% Monodromy matrix
Mo = reshape(Y(end, 7:42), 6, 6);

% Eigenvector and eigenvalue analysis
[V, D] = eig(Mo); %dim=6*6
disp(diag(D))
for i=1:d
    if imag(D(i,i))>1e-3
        idxcen=i;
    end
end

% Eigenvalue corresponding to center direction
Ec = D(idxcen,idxcen);

% Eigenvector corresponding to center direction
Vc = V(:, idxcen);
if Vc(1)<0
    Vc = -Vc;
end

% rotation number
% (Ref:(p162))
rho = atan2(imag(Ec),real(Ec));

% replicate and reshape the vector containing winding angles
THTj = repmat(thtj,d,1);
THTj = reshape(THTj,[],1);

% perturbation
x_star = zeros(d,M);
uhat = zeros(d*M,N);

% Correction factor
%(Ref:(5))
cf = exp(-1i*rho*tk/Tpo);

for idx=1:M
    % vec index
    vdx = d*idx-(d-1):d*idx;

    % time at index
    t = tk(idx);

    % STM
    if idx==1
        x_star(:,idx) = xpo;
        phit = eye(d);
    else
        [~,Y] = ode113(@(t,x) fun_stm_cr3bp(t,x,mu),[0.0 t],XP0,options);
        x_star(:,idx) = Y(end,1:d)';
        phit = reshape(Y(end,d+1:d+d^2),d,d);
    end

    % propagate
    %(Ref:(5))
    vi = cf(idx)*phit*Vc;
    VI = repmat(vi,N,1);

    % perturbation
    %(Ref:(4))
    uh = K.*(cos(THTj).*real(VI)-sin(THTj).*imag(VI));
    uhat(vdx,:) = reshape(uh,d,N);
end
x_star = repmat(reshape(x_star,[],1),1,N);
x_star = x_star + uhat;%[d*M,N]

% organize curves
if p("ctswitch")
    x_starg = zeros(d*N*M,1);
    uhatg = zeros(d*N*M,1);
    for i=1:M
        k=d*N;
        idx = i*k-(k-1):i*k;
        % invariant curve
        u = x_star(d*i-(d-1):d*i,:);
        % store invariant curve
        x_starg(idx) = u;
        % tangent direction
        uh = uhat(d*i-(d-1):d*i,:);
        % store tangent direction
        uhatg(idx) = uh;
    end
    U = x_starg;
    uhat = uhatg;
end

% frequencies
w = [2*pi/Tpo;rho/Tpo];

% fourier matrices
DFr = FR{3,2};

% initial angle partials
DU0t1 = DFr*U(1:d*N,1);
DU0t0 = zeros(d*N,1);
for j=1:N
    f = fun_cr3bp([],U(d*j-(d-1):d*j),mu);
    DU0t0(d*j-(d-1):d*j,1) = (1/w(1)).*(f-w(2).*DU0t1(d*j-(d-1):d*j,1));
end

% QPO initial guess
z0 = [U;Tpo;rho];
% tangent vector
phi0 = [uhat;zeros(2,1)]/sqrt(dot(uhat,uhat)/N*M);

% QPO solution data package
Ud0 = cell(7,2);
% initial angle partials
Ud0{1,1} = "DU0t0";
Ud0{1,2} = DU0t0;
Ud0{2,1} = "DU0t1";
Ud0{2,2} = DU0t1;
% initial solution
Ud0{3,1} = "z0";
Ud0{3,2} = z0;
% patch times
Ud0{4,1} = "pt";
Ud0{4,2} = linspace(0,Tpo-1e-5,M+1);
Ud0{5,1} = "Ut";
Ud0{5,2} = zeros(d*N,M);
% initialize matrix of STM
Ud0{6,1} = "phiT";
Ud0{6,2} = zeros(d*N,d*N,M);
% vector field at final integrated state
Ud0{7,1} = "fT";
Ud0{7,2} = zeros(d*N,1);

end