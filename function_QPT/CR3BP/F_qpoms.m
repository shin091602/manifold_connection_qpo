function [F,Ud] = F_qpoms(Z,zpo,Ud,p,FR)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute the error vector in calculating QPO by using GMOS method
%% By:Soi Yamaguchi

%%% input
% Z :current state variables [X;T;rho;w]
% Ud :torus solution package
% p :parameter dictionary
% FR :fourier matrices

%%% output
% F :error vector
% Ud :updated torus solution package
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% DICTIONARY OPEN
d = p("d");
N = p("N");
M = p("M");
mu = p("mu");

%% OPTIONS ODE
options = odeset('RelTol',1e-13, 'AbsTol',1e-13);

% current state variables
X = Z(1:d*N*M);
T = Z(d*N*M+1);
rho = Z(d*N*M+2);

% fourier matrices
Fr = FR{1,2};
IFr = FR{2,2};
Q = FR{4,2};
Q = rotation_matrix(Q,rho,p);
R = IFr*Q*Fr;

%% INTEGRATION
% initialization matrix of integrated states
Ud{5,2} = zeros(d*N,M);
% initialize matrix of STM
Ud{6,2} = zeros(d*N,d*N,M);
% integrated time
inti = [Ud{4,2}';T];
for i=1:M
    % designate invariant circle
    k=d*N;
    U = X(i*k-(k-1):i*k);
    % initialize final state
    Uf = zeros(d*N,1);
    for j=1:N
        % designate initial state
        u = U(d*j-(d-1):d*j);
        % time span
        ts = [inti(i) inti(i+1)];
        % ODE
        [~,Y] = ode113(@(t,x) fun_stm_cr3bp(t,x,mu),ts,[u;reshape(eye(d),[],1)],options);
        % final state
        Uf(d*j-(d-1):d*j) = Y(end,1:d)';
        % monodromy matrix
        phiT = reshape(Y(end,d+1:d+d^2),d,d);
        % state STM
        Ud{6,2}(d*j-(d-1):d*j,d*j-(d-1):d*j,i) = phiT;
    end
    % integrated states
    Ud{5,2}(:,i) = Uf;
end
% vector field at final state
UdT = Ud{5,2}(:,end);
fT = zeros(d*N,1);
for j=1:N
    fT(d*j-(d-1):d*j) = fun_cr3bp([],UdT(d*j-(d-1):d*j),mu);
end
Ud{7,2} = fT;

%define error vector
F = zeros(d*N*M+p("pcon"),1);

%% INVARIANCE CONDITION
%(Ref:(21))
F(1:d*N) = R*UdT-X(1:d*N);

%% MULTIPLE SHOOTING CONDITIONS
for i=2:M
    %invariant circle index
    k=d*N;
    %patch points
    uptch = X(i*k-(k-1):i*k);
    %difference between integrated previous points
    F(i*k-(k-1):i*k) = Ud{5,2}(:,i-1)-uptch;
end

%% PHASE CONDITIONS
%(Ref:(p171 p_0 and p_1))
% previous solution
zpre = Ud{3,2};
% initial angle partials
DU0t0 = Ud{1,2};
DU0t1 = Ud{2,2};
% calculate phase condition
p0=(1/N).*(X(1:d*N)-zpre(1:d*N))'*DU0t0;
p1=(1/N).*(X(1:d*N)-zpre(1:d*N))'*DU0t1;
F(d*N*M+1:d*N*M+2) = [p0;p1];

%% FIX STROBOSCOPIC TIME
%(Ref:(22))
s0z = Z(d*N*M+1)-zpo(end);

%% ERROR VECTOR
F = [F;s0z];

end