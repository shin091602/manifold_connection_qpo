% CR3BP
% matrix ver
% For 2D-tori
% For plot manifold
% For plot torus

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% calculate 2D quasi-periodic invariant tori (QPT) by leveraging GMOS algorithm
%% matrix version
%% created by Soi Yamaguchi & Ayano Tsuruta & Soichiro Shin
%% updated in 4/21/2025
%% based on Julia written by Damennick Bolte Henry

%% Outputs:
%sol_qpos:converged solution packages

%% Textbook:
%Nicola Baresi et al., "Fully Numerical Methods for Continuing Families of Quasi-Periodic Invariant Tori in Astrodynamics"
%Refer to pages:p169-p172
%Notation rules:
%(Ref:(num)) -- refer to the equation number in the textboook.

%% Explanation:
%Refer to public5/M2(2023年度)/山口/workspace/GMOS_algorithm.pptx

%% Functions:
%po_mulshoot_initialization.m
%F_poms.m
%DF_poms.m
%PAC_poms.m
%PAC_correction_po.m
%Xd_finalization_poms.m
%plot_refined_orbit.m
%plot_invariant_curve.m
%PAC_qpoms.m
%PAC_correction_qpo.m
%F_GMOS.m
%DF_GMOS.m
%rotation_matrix.m
%fourier_matrix.m
%Xd_finalization_qpoms.m
%plot_qpos.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;  %close all figures
clear;      %clear all variables
clc;        %clear the command terminal
format long
%warning off

% line width
set(0, 'DefaultLineLineWidth', 1.2) % default 0.5pt
set(0, 'DefaultAxesLineWidth', 1.2)
set(0, 'DefaultTextLineWidth', 1.2)

% font size
set(0, 'DefaultTextFontSize', 24)
set(0, 'DefaultAxesFontSize', 24)

% font name
set(0, 'DefaultTextFontName', 'Times New Roman')
set(0, 'DefaultAxesFontName', 'Times New Roman')
set(0, 'DefaultTextInterpreter', 'Latex')
set(0, 'DefaultLegendInterpreter', 'Latex')

% figure color
set(0, 'DefaultFigureWindowStyle', 'docked');
set(gcf, 'Color', 'none');
set(gca, 'Color', 'none');
set(gcf, 'InvertHardCopy', 'off');

close

%% addpath
current_pass = pwd;

% BASIC
addpath(append(current_pass, '/Ayano functions'));
addpath(append(current_pass, '/Ayano functions/A'));
addpath(append(current_pass, '/Ayano functions/ODE/natural'));
addpath(append(current_pass, '/Ayano functions/ODE/nonnatural'));
addpath(append(current_pass, '/Ayano functions/plot'));
addpath(append(current_pass, '/Functions'));
% GMOS
addpath(append(current_pass, '/function_QPT'));
addpath(append(current_pass, '/function_QPT/CR3BP'));
% SUSUMU
addpath(append(current_pass, '/susumu functions'));
myTimer = tic;        %start timer

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% INITIAL VALUES
% L2 halo orbit --a sample orbit
x0 = [1.116913489266650; 0; 0.021776489613409; 0; 0.186114812512697; 0];%state
T = 3.407900091010860;%period

%% DICTIONARY OF THE VARIABLES
%%% notes:
% constant values for the calculation are defined here.
p = dictionary();

p("Ms") = 5.9724e+24;%mass of Sun　→ earth
p("M1") = 7.3458e+22;%mass of Earth　→ moon
p("mu") = 1.21536E-02;%mu--EM
p("G") = 6.67300000000000e-11;%gravitational constant
p("chara_length_CR3BP") = 3.844e+8;%[m]
p("chara_mass_CR3BP") = 5.9724e+24;%[kg]
p("chara_time_CR3BP") = 3.7748e+5;%[s]
p("N_CR3BP_SE") = 1.6645e-05;

p("d") = 6; % dimension of variables
p("N") = 15; % number of cross section of invariant circle. 3とかだとガタガタになる．
p("M") = 15; %number of invariant circle --the number of N has to be equal to the number of M.→離散フーリエ変換の都合上．あとM,Nは奇数じゃないとだめ．

p("K") = 2e-2; %%radius of invariant circle.大きすぎるとcontinuationできない．小さすぎると，ランクが下がる．initiral = 1e-4

p("Iteration") = 50; %iteration of GMOS　%初期設定は50 % continuationの回数1000回にしてもfamilyほとんど変わらないが，少し高さが違う．．
p("Threshold") = 1e-8; %converge threshold (初期設定は1e-8)
p("C") = fun_Jacobi_const_CR3BP(x0, p("mu")); %Jacobi constant
p("phsflg") = true;%flag of phase condition --ON
p("ctswitch") = true;%center tangenet switch
p("pcon") = 2;%phase condition index
p("congflg") = false;%converged flag
p("HLflg") = true; %Halo:true, Lyapunov:false
p("rho_type") = false; %true:in-plane, false:out-plane

p("pert") = 1e-4;% perturbation
%** for manifold **%
p("snap_ini_time") = 0.8*pi;%initial time of snapshots %default:pi/2
p("snap_fin_time") = 1.2*pi;%final time of snapshots %default:0.8pi
% p("snap_ini_time") = 0.5*pi;%initial time of snapshots %default:pi/2
% p("snap_fin_time") = 2.5*pi;
p("snap_span") = 2;%number of snapshots defalut:6 (snapを何個取るか．多分manifoldのとき)
p("num_iter") = 1080;%number of grid points for interpolation initial: 200
snap_time = linspace(p("snap_ini_time"),p("snap_fin_time"),p("snap_span"));

%% OPTIONS ODE
options_ODE = odeset('RelTol',1e-13, 'AbsTol',1e-13, 'Refine', 10);

%% LIBRATION POINTS
[~,L2,~,~,~] = librationPoints(p("mu"));

% PO
[t_PO, x_PO] = ode113(@(t,x) fun_ode_n_CR3BP(t, x, p("mu")),[0 T],x0,options_ODE);
x_PO = x_PO';

%% plot modromy's eigenvalue (in-plane:blue, out-plane:red, in+out:green)
% Monodromy matrix
% initial state including initial STM
XP0 = [x0; reshape(eye(p("d")), [], 1)];
[~,Y] = ode113(@(t,x) fun_stm_cr3bp(t,x,p("mu")),[0 T],XP0,options_ODE);

% Monodromy matrix
Mo = reshape(Y(end, p("d")+1:p("d")*p("d")+p("d")), p("d"), p("d"));

% Eigenvector and eigenvalue analysis
[V_n, D_n] = eig(Mo); %dim=6*6

[~, fig_monodromy_eig_zoom] = funplot_monodromy_eigenvalue(diag(D_n), V_n);

%% PATCH POINTS -- center of invariant curves（断面上の青い点を決めている）
z0 = [x0;T];
[zpo,Xd] = po_mulshoot_initialization(z0,p); % 断面ごとの周期軌道の座標[r;v](zpo(1:6*M))と周期(zpo(end))を返す．Xdの意味は単なる確認用．
[hPO,orbitPO] = plot_periodic_orbit(zpo,p,L2);

%% CALCULATE INITIAL INVARIANT CIRCLE(=周期軌道の断面) AND ROTATION NUMBER（monodromyのeigenvalueを使って，ρの初期値を出している．線形な範囲では，固有値が1ではない周期解のarctan(Im/Re)がρ）．意味としては，"w0に対し，どれだけ振動するか，位相がずれているか"ｇρの定義．非線形の範囲までcontinuationしても閉じるようにdifferential correctionが行われる．
[z0,phi0,Ud0,rho,Ec,x_star_pert] = fun_center_manifold_CR3BP_matrix([zpo(1:p("d"));zpo(end)],p);
% z0が変数[V0, ρ, T]の初期推定解(M個の円の，ｎ分割分の解)
[hinvcir, check_accuracy] = plot_invariant_curve_CR3BP_matrix(x_star_pert, zpo,z0,p); % ちゃんと-ρしたらもとに戻ったか
C_periodic = fun_Jacobi_const_CR3BP(zpo(1:6), p("mu")); % Jacobi constant of periodic orbit, almost same with C_periodic_original

%% DICTIONARY OF PSEUDO-ARCLENGTH CONTINUATION
%%% notes:
% additional constant values for the pseudo-arclength continuation are defined here.
pacqp = dictionary();
pacqp("n") = 5;%iteration of PAC --if one need to calculate the family, change its value. %初期設定5
pacqp("tol") = 5e-7;%error tolerance
pacqp("itmax") = 100;%max iteration
pacqp("optit") = 5;%maximum step length
pacqp("smax") = 1e-5;%maximum step length
pacqp("jcr") = p("d")*p("N")*p("M")+p("pcon");%Jacobian rows to consider for initial family tangent
pacqp("issparse") = true;%sparse Jacobian flag
pacqp("N") = p("N")*p("M");%number of continuation function points
pacqp("fpidx") = p("d")*p("N")*p("M");%continuation function point indices
pacqp("plotflg") = false;%plot flag
pacqp("coflg") = false;%correction only flag
pacqp("fdcheck") = false;%finite difference check

%% PSEUDO-ARCLENGTH CONTINUATION
% fix T: sol_qpos = PAC_qpoms_CR3BP_matrix(z0,zpo,1e-5,Ud0,phi0,p,pacqp);
sol_qpos = PAC_qpoms_CR3BP_matrix(z0,zpo,1e-5,Ud0,phi0,p,pacqp);%(Ref:(19)-(23))

% fix Jacobi constant
% sol_qpos = PAC_qpoms_CR3BP_matrix_Jacobi_constant_fix(z0,zpo,1e-5,Ud0,phi0,p,pacqp,C_periodic);%(Ref:(19)-(23))

%% save results
% data_name = strcat('Natural_Lyapunov_2D_matrix_inplane_QPT_CR3BP_L2_EM');
% data_name = strrep(data_name,'.',',');
% save(data_name);

%% Focus on the last torus and interpolate
%solution
solc_qpos = sol_qpos{pacqp("n"),1};
fin_qpos = solc_qpos{1,2};

% INTERPOLATION（2025点分のデータを40000点に拡張）
% collect grid points
% assuming N==M, meshgrid
%tht0, tht1 :meshgrid
tht0 = linspace(-2*pi,4*pi,p("N")*3); % dim=45
tht1 = linspace(-2*pi,4*pi,p("M")*3)'; % dim=45
%tht0_n,tht1_n :interpolated grid
tht0_n = linspace(-2*pi,4*pi,p("num_iter")); % dim=200
tht1_n = linspace(-2*pi,4*pi,p("num_iter")); % dim=200
% interpolate
[ri,~] = interpolate_torus(fin_qpos,tht0,tht1,tht0_n,tht1_n,p); %実データfin_qpos(dim=225)から，interpolate_torusも用いて間を埋める（225点のままでは粗く，曲面を滑らかにするため．）．40000点に拡張する．[th0, th1]→ri=[X,Y,Z,Xdot,Ydot,Zdot]. Xは200*200の行列．XはX_1からX_40000まである．


%% FOCUS ON ONE TORUS
% designate a torus
nm = 5;%1~pacqp("n")
solc_qpos_m = sol_qpos{nm,1};
fin_qpos_m = solc_qpos_m{1,2};

%% MANIFOLD SURF PLOT
% calculate the directions of the designated manifolds (-ρしたmonodromyの固有ベクトルを使ってmanifoldのdirection決める)
del_w_us_inter_full = direction_manifold_CR3BP_matrix(fin_qpos_m,p);

% calculate the manifolds orientation of the interpolated points
del_w_us_interpolation = direction_manifold_interpolation(fin_qpos_m,p);

% vector arrows
U0 = zeros(p("d")*p("M"),1);

% invariant curve idx
for i=1:p("M") %designate the invariant circle
    k=p("d")*p("N")*(i-1)+1:p("d")*p("N")*i;
    X0 = fin_qpos_m(k);
    del_w_us = del_w_us_inter_full(:,:,i);
    for j=1:p("N")%designate torus point
        % Initial points of torus unstable manifolds
        U0(p("d")*i-(p("d")-1):p("d")*i) = X0(p("d")*j-(p("d")-1):p("d")*j)+p("pert").*del_w_us(:,j);
    end
end

% plot
hinterm = figure();
hold on
grid on
box on
axis equal
view([7 25]);
xlabel('$x$[-]');
ylabel('$y$[-]');
zlabel('$z$[-]');


%calculate ZVC
[x,y] = meshgrid(-1.5:1e-3:1.5);
[a,b] = size(x);
z = zeros(a,b);

x = reshape(x, [a*b,1]);
y = reshape(y, [a*b,1]);
z = reshape(z, [a*b,1]);

r = sqrt(x.^2+y.^2+z.^2);
U = 1./r + 1/2*(3*x.^2-z.^2);
c = 2.*U;

x = reshape(x, [a,b]);
y = reshape(y, [a,b]);
z = reshape(z, [a,b]);
c = reshape(c, [a,b]);

% surf torus
hs = surf(ri{1,1},ri{2,1},ri{3,1},CO);
hold on

% sample trajectories of torus manifolds
% 1th and 2nd points are representives
for i=1:5:p("N")
    [~,X_sample] = ode113(@(t,x) fun_ode_n_CR3BP(t, x, p("mu")),[0 p("snap_fin_time")+0.1],U0(p("d")*i-(p("d")-1):p("d")*i));
    plot3(X_sample(:,1),X_sample(:,2),X_sample(:,3),"r","LineWidth",1.5);
    hold on
end

% interpolate function
[rim,~,~,~] = interpolate_manifold(fin_qpos_m,del_w_us_inter_full,tht0,tht1,tht0_n,tht1_n,p);

% surf manifold
%color
CO_un = zeros(p("num_iter"),p("num_iter"),3);
CO_un(:,:,1) = 0.7010.*ones(p("num_iter")); % red

% initialization of frame
for sn=1:p("snap_span")
    surf(real(rim{1,sn}),real(rim{2,sn}),real(rim{3,sn}),CO_un);
end

% for visualization
shading interp
lightangle(27,36)
lightangle(27,36)
hold off

%% End of script
time = strcat('calculation time: ', num2str(toc(myTimer)));
disp(time);

