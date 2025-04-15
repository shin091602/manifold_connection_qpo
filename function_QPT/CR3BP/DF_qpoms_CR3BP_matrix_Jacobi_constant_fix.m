function [DF,Ud] = DF_qpoms_CR3BP_matrix_Jacobi_constant_fix(Z,Ud,p,zpo)
% function [DF,Ud] = DF_qpoms_2n(Z,Ud,p, zpo,R_eye, Qup, Qdown, Q_1313, Q_1346, Ja, xt, F)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute the error Jacobian in calculating QPO by using GMOS method
%%% input
% Z :current state variables [X;T;rho;w]
% Ud :torus solution dictionary

%% QPO solution data package
% Ud0{1,1} = "DU0t0"; % トーラス関数uを角度theta_0でそれぞれ偏微分したもの。
% % % 目的関数内にあるp0とかp1を求めるときに使うんだけど，correctionする前に最初のiterationでの目的関数Fの値を求める必要があるから，fun_center_manifold~の関数内であらかじめ計算してる
% Ud0{1,2} = DU0t0;
% Ud0{2,1} = "DU0t1"; % トーラス関数uを角度theta_1でそれぞれ偏微分したもの。
% Ud0{2,2} = DU0t1;
% % initial solution
% Ud0{3,1} = "z0";
% Ud0{3,2} = z0;
% % patch times
% Ud0{4,1} = "pt";
% Ud0{4,2} = linspace(0,Tpo-1e-5,M+1);
% Ud0{5,1} = "Ut"; = 1T後の座標yf
% Ud0{5,2} = zeros(d*N,M);
% % initialize matrix of STM
% Ud0{6,1} = "phiT";
% Ud0{6,2} = zeros(d*N,d*N,M);
% % vector field at final integrated state
% Ud0{7,1} = "fT"(M=1は先頭に，M=15は最後に格納している);
% Ud0{7,2} = zeros(d*N,1);

% p :parameter dictionary
% FR :fourier matrices

%%% output
% F :error vector
% Ud :updated torus solution package
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% DICTIONARY OPEN
d = p("d");
M = p("M");
N = p("N");

% current state variables
rho = Z(d*N*M+2);

% % fourier matrices
[FR, ~, DFR_rho] = fun_Fourier_matrix(rho,p);

% initialize Jacobian
DF = zeros(d*N*M+p("pcon"),d*N*M+p("pcon"));

%% INVARIANCE CONDITION PARTIALS
% wrt U0
DF(1:d*N,1:d*N) = -eye(d*N,d*N);
DF(1:d*N,end-(d*N)-1:end-2) = FR*Ud{6,2}(:,:,end);

% wrt T
DF(1:d*N,d*N*M+1) = FR*Ud{7,2}(1:d*N,1);

% wrt rho
DF(1:d*N,d*N*M+2) = DFR_rho*Ud{5,2}(:,end);

%% MULTIPLE SHOOTING PARTIALS
for i=2:M
    k=d*N;
    % row idx
    ridx = (i*k-(k-1)):i*k;
    % column idx
    cidx = (i*k-(k-1))-k:i*k-k;
    % partials
    DF(ridx,cidx) = Ud{6,2}(:,:,i-1);
    DF(ridx,ridx) = -eye(k);
    DF(ridx,d*N*M+1) = Ud{7,2}(cidx,1).*(i-1)/M; %必要な気がするけど，，，
end

%% PHASE CONDITION PARTIALS
DP = [(1/N).*Ud{1,2}';(1/N).*Ud{2,2}'];
DF(d*N*M+1:d*N*M+2,1:d*N) = DP(1:2,:);

%% dC
for i = 1:M
    for j = 1:N
        vec_begin = d*N*(i-1)+d*(j-1)+1;
        vec_end = d*N*(i-1)+d*j;
        dC(vec_begin:vec_end, 1)= Jacobi_const_diff_CR3BP(Z(vec_begin:vec_end, 1),p("mu"));
    end
end
Dsi = [dC',0,0];

DF = [DF;Dsi];

% %% check the Finite difference
% Z_check = cell(1,d*M*N+4);
% Z_check{1} = Z + [1e-8; zeros(d*N*M+3,1)];
% Z_check{d*N*M+4} = Z + [zeros(d*N*M+3,1); 1e-8];
% 
% for i = 2:d*M*N+3
%     Z_check{i} = Z + [zeros(i-1,1); 1e-8; zeros(d*N*M+4-i,1)];
% end
% 
% F_check = cell(1,d*M*N+4);
% check_finite = cell(1,d*M*N+4);
% check_finite_norm = zeros(1,d*M*N+4);
% 
% [F,~] = F_qpoms_Hill3BP_2n(Z,zpo,Ud,p, R_eye, Qup, Qdown, Q_1313, Q_1346, Ja, xt);
% for i=1:d*N1*N2*2
%     [F_check{i},~,~] = F_qpoms_Hill3BP_2n_3D_newnew(Z_check{i},zpo,Ud,p, R_eye, Qup, Qdown, Q_1313, Q_1346, Ja, xt);
%     a = (F_check{i}-F)/1e-8;
%     check_finite{i} = DF(:,i)-a;
%     check_finite_norm(i) = norm(check_finite{i});
% end
% 
% for i=1:d*N*2
%     [F_check{i},~] = F_qpoms_Hill3BP_2n_newnew(Z_check{i},zpo,Ud,p, R_eye, Qup, Qdown, Q_1313, Q_1346, Ja, xt);
%     a = (F_check{i}-F)/1e-8;
%     check_finite{i} = DF(:,i)-a;
%     check_finite_norm(i) = norm(check_finite{i});
% end
% 
% 
% for i=d*N*M+1:d*N*M+4
%     [F_check{i},~] = F_qpoms_Hill3BP_2n_newnew(Z_check{i},zpo,Ud,p, R_eye, Qup, Qdown, Q_1313, Q_1346, Ja, xt);
%     a = (F_check{i}-F)/1e-8;
%     check_finite{i} = DF(:,i)-a;
%     check_finite_norm(i) = norm(check_finite{i});
% end
end