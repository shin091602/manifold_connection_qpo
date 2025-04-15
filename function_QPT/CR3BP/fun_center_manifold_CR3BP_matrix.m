function [z0,phi0,Ud0,rho,Ec, x_star_pert] = fun_center_manifold_CR3BP_matrix(zpo,p)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% calculate invariant circle and rotation number

%%% input
%zpo :initial data of period
%p :parameter dictionary
%FR :fourier matrix dictionary

%%% output
%z0 :initial guess of invariant circle
%phi0 :initial tangent vector
%Ud0 :torus function dictionary

%% flow
% まず，周期軌道を初期値から1周期伝播して，そのモノドロミー行列の固有値を使って一つ目の不変円を作る．
% 次に，1周期をM分割したときの各時間における状態量をそこでの不変円の作成に使ってた.
% だから周期軌道の設計にはmultiple shootingは利用してない！
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% DICTIONARY
d = p("d");
M = p("M");
N = p("N");
K = p("K");
mu = p("mu");

%% OPTIONS ODE
options_ODE = odeset('RelTol',1e-13, 'AbsTol',1e-13);

%% calculate eigenvalues of the monodromy matrix
% Propogate orbit for one period and gather data at N points
% periodidx_icen orbit state
xpo = zpo(1:d); % 1個目の断面の初期値しか持ってきてない．
% period
Tpo = zpo(end);
% time vector % (dim=1*15)
t_i = linspace(0,Tpo,M+1); % dim=1*15
t_i = t_i(1,1:end-1);
% invariant circle angles (dim=1*15)
thtj = linspace(0,2*pi,N+1); %invariant circleをN分割
thtj = thtj(1,1:end-1);

% initial state including initial STM
XP0 = [xpo; reshape(eye(d), [], 1)]; 
[~,Y] = ode113(@(t,x) fun_stm_cr3bp(t,x,mu),[0 Tpo],XP0,options_ODE);

% Monodromy matrix
Mo = reshape(Y(end, d+1:d*d+d), d, d);

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
    Vc = -Vc; % xを正の向きのベクトルにする, やってもやらなくてもいい
end

% rotation number
% (Ref:(p162))
rho = atan2(imag(Ec),real(Ec));

% % perturbation
x_star = zeros(d,M);
% uhat = zeros(d*M,N);

% Correction factor
%(Ref:(5))
cf = exp(-1i*rho*t_i/Tpo); % dim=1*15

phij = zeros(d,N,M);
phij_i = zeros(d,N,M);
x_star_pert = zeros(d,N,M);
% ここからM分割する
for idx_i=1:M
    % STM
    if idx_i==1
        x_star(:,idx_i) = xpo;
        phit = eye(d);
    else
        [~,Y] = ode113(@(t,x) fun_stm_cr3bp(t,x,mu),[0.0 t_i(idx_i)],XP0,options_ODE);
        x_star(:,idx_i) = Y(end,1:d)';
        phit = reshape(Y(end,d+1:d+d^2),d,d);
    end
    %% とりあえず2D verをmatrix表記に変換
    for idx_j = 1 : N
        % invariant circle, eq(3.3)
        phij(:, idx_j, idx_i) = K.*(cos(thtj(idx_j)).*real(Vc*cf(idx_i))-sin(thtj(idx_j)).*imag(Vc*cf(idx_i))); %2Dver : dim(phij) = d*N
        % perturbation % eq(3.4) %(Ref:(5))
        phij_i(:, idx_j, idx_i) = phit*phij(:, idx_j, idx_i); % =uh, dim = 12*N*M
        x_star_pert(:, idx_j, idx_i) = x_star(:, idx_i) + phij_i(:, idx_j, idx_i);
        % x_star_pert_tilda{i}(j1,j2,:) = x_star_pert{j1, j2, i};
    end
end
% x_star = repmat(reshape(x_star,[],1),1,N);
% x_star = x_star + uhat; %[d*M,N]=[90*15]

% organize curves　1からMについてそれぞれN個つくる
if p("ctswitch")
    % x_starg = zeros(d*N*M,1);
    % uhatg = zeros(d*N*M,1);
    % for i=1:M
    %     k=d*N;
    %     idx_i = i*k-(k-1):i*k; %dim=1*90=6*NでiはM
    %     % invariant curve (M個の不変円)
    %     u = x_star(d*i-(d-1):d*i,:); % dim=6*15=6*N
    %     % store invariant curve
    %     x_starg(idx_i) = u; %dim=6*15(N)を15個(M)繋げる
    %     % tangent direction (ρの方向)　dim=6*15
    %     uh = uhat(d*i-(d-1):d*i,:);
    %     % store tangent direction
    %     uhatg(idx_i) = uh;
    % end
    % U = x_starg; % dim=1350*1
    % uhat = uhatg; % dim=1350*1
    U = reshape(x_star_pert, d*N*M, 1);
    uhat = reshape(phij_i, d*N*M, 1);
end

%% fourier matrix
[~, DFR_theta, ~] = fun_Fourier_matrix(rho,p);

% frequencies
w = [2*pi/Tpo; rho/Tpo];

%% For Phase constraint
% for g = 1:d
%     for j = 1:N
%         DU0t1(:,j) = DDFT(:,:,g)*x_star_pert(g,j,1); % dim(DUOt1)=90 %x_star_pert = zeros(d,N,M);
%     end
% end

% initial angle partials
DU0t1 = DFR_theta*U(1:d*N,1); % dim(DUOt1)=180
DU0t0 = zeros(d*N,1); % dim(DUOt0)=180
for j=1:N
    f = fun_cr3bp([],U(d*j-(d-1):d*j),mu);
    DU0t0(d*j-(d-1):d*j,1) = (1/w(1)).*(f-w(2).*DU0t1(d*j-(d-1):d*j,1));
end

% QPO initial guess
z0 = [U;Tpo;rho]; % dim(z0)=1352(U=1350, dim(T),dim(rho)=1)
% tangent vector
phi0 = [uhat;zeros(2,1)]/sqrt(dot(uhat,uhat)/N*M); % dim(phi0)=1352

% QPO solution data package
Ud0 = cell(7,2);
% % initial angle partials
Ud0{1,1} = "DU0t0"; % トーラス関数uを角度theta_0でそれぞれ偏微分したもの。
% % 目的関数内にあるp0とかp1を求めるときに使うんだけど，correctionする前に最初のiterationでの目的関数Fの値を求める必要があるから，fun_center_manifold~の関数内であらかじめ計算してる
Ud0{1,2} = DU0t0;
Ud0{2,1} = "DU0t1"; % トーラス関数uを角度theta_1でそれぞれ偏微分したもの。
Ud0{2,2} = DU0t1;
% initial solution
Ud0{3,1} = "z0";
Ud0{3,2} = z0;
% patch times
Ud0{4,1} = "pt";
Ud0{4,2} = linspace(0,Tpo,M+1);
%Ud0{4,2} = linspace(0,Tpo-1e-5,M+1);
Ud0{5,1} = "Ut";
Ud0{5,2} = zeros(d*N,M);
% initialize matrix of STM
Ud0{6,1} = "phiT";
Ud0{6,2} = zeros(d*N,d*N,M);
% vector field at final integrated state
Ud0{7,1} = "fT";
Ud0{7,2} = zeros(d*N,1);

end