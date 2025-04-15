% CR3BP
% 2D-tori

% 2024.10.24
% rhoにマイナスはつけてない
% wは変数からはずした

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% calculate 2D quasi-periodic invariant tori (QPT) by leveraging GMOS algorithm
%% modified by Ayano Tsuruta (2024)
%% created by Soi Yamaguchi (2023)
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
addpath(append(current_pass, '\Ayano functions'));
addpath(append(current_pass, '\Ayano functions\A'));
addpath(append(current_pass, '\Ayano functions\ODE\natural'));
addpath(append(current_pass, '\Ayano functions\ODE\nonnatural'));
addpath(append(current_pass, '\Ayano functions\plot'));
addpath(append(current_pass, '\Functions'));
% GMOS
addpath(append(current_pass, '\function_QPT'));
addpath(append(current_pass, '\function_QPT\Hill3BP'));
addpath(append(current_pass, '\function_QPT\CR3BP'));
myTimer = tic;        %start timer

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% INITIAL VALUES
% L2 halo orbit --a sample orbit
x0 = [1.010059516040882;0;0.008126717418058;0;-0.013983743043142;0];%state
T = 3.0031;%period

%% DICTIONARY OF THE VARIABLES
%%% notes:
% constant values for the calculation are defined here.
p = dictionary();
p("Ms") = 1.9885e+30;%mass of Sun
p("M1") = 5.9724e+24;%mass of Earth
p("mu") = 3.00346093142064e-06;%mu--SE
p("G") = 6.67300000000000e-11;%gravitational constant
p("chara_length_CR3BP") = 149600000000.000;%[m]
p("chara_mass_CR3BP") = 1.98850597240000e+30;%[kg]
p("chara_time_CR3BP") = 5023117.97170552;%[s]
p("N_CR3BP_SE") = 1.99079536979392e-07;
p("d") = 6; % dimension of variables
p("N") = 15; % number of points in circle
p("M") = 15; %number of circle --the number of N has to be equal to the number of M.
p("K") = 2e-3; %radius of circle
p("null_acr") = p("K")*1e-1; %null accuracy
p("Iteration") = 50; %iteration of GMOS
p("Threshold") = 1e-8; %converge threshold
p("C") = Jacobi_const(x0,p("mu"));%Jacobi constant
p("phsflg") = true;%flag of phase condition --ON
p("ctswitch") = true;%center tangenet switch
p("pcon") = 2;%phase condition index
p("congflg") = false;%converged flag
p("HLflg") = true;%Halo:true, Lyapunov:false
p("pert") = 1e-4;%deorbit velocity
%** for manifold **%
p("snap_ini_time") = 0.8*pi;%initial time of snapshots %default:pi/2
p("snap_fin_time") = 1.2*pi;%final time of snapshots %default:0.8pi
p("snap_span") = 2;%number of snapshots defalut:6
p("num_iter") = 200;%number of grid points for interpolation
snap_time = linspace(p("snap_ini_time"),p("snap_fin_time"),p("snap_span"));

%% OPTIONS ODE
options_ODE = odeset('RelTol',1e-13, 'AbsTol',1e-13);

%% LIBRATION POINTS
[~,L2,~,~,~] = librationPoints(p("mu"));

%% PATCH POINTS -- center of invariant curves
z0 = [x0;T];
[zpo,Xd] = po_mulshoot_initialization(z0,p);
[hPO,orbitPO] = plot_periodic_orbit(zpo,p,L2);

%% FOURIER MATRICES PACKAGE
[R,IR,DR] = fourier_matrix(p);
FR = cell(5,2);
FR{1,1} = "Fr";
FR{1,2} = R;
FR{2,1} = "IFr";
FR{2,2} = IR;
FR{3,1} = "DFr";
FR{3,2} = DR;
FR{4,1} = "Q";
FR{4,2} = zeros(p("d")*p("N"),p("d")*p("N"));
FR{5,1} = "R";
FR{5,2} = zeros(p("d")*p("N"),p("d")*p("N"));

%% CALCULATE INITIAL INVARIANT CIRCLE AND ROTATION NUMBER
[z0,phi0,Ud0,rho,Ec] = fun_center_manifold_cr3bp([zpo(1:p("d"));zpo(end)],p,FR);%(Ref:(4)-(5))
[hinvcir] = plot_invariant_curve(zpo,z0,p,FR);

%% DICTIONARY OF PSEUDO-ARCLENGTH CONTINUATION
%%% notes:
% additional constant values for the pseudo-arclength continuation are defined here.
pacqp = dictionary();
pacqp("n") = 5;%iteration of PAC --if one need to calculate the family, change its value.
pacqp("tol") = 1e-7;%error tolerance
pacqp("itmax") = 100;%max iteration
pacqp("optit") = 5;%maximum step length
pacqp("smax") = 1e-5;%maximum step length
pacqp("jcr") = p("d")*p("N")*p("M")+4;%Jacobian rows to consider for initial family tangent
pacqp("issparse") = true;%sparse Jacobian flag
pacqp("N") = p("N")*p("M");%number of continuation function points
pacqp("fpidx") = p("d")*p("N")*p("M");%continuation function point indices
pacqp("plotflg") = false;%plot flag
pacqp("coflg") = false;%correction only flag
pacqp("fdcheck") = false;%finite difference check

%% PSEUDO-ARCLENGTH CONTINUATION
sol_qpos = PAC_qpoms(z0,zpo,1e-5,Ud0,phi0,p,pacqp,FR);%(Ref:(19)-(23))

%% TORUS FAMILY
for n=1:pacqp("n")
  %solution
  solc_qpos = sol_qpos{n,1};
  fin_qpos = solc_qpos{1,2};

  %% INTERPOLATION
  % collect grid points
  % assuming N==M, meshgrid
  %tht0, tht1 :meshgrid
  tht0 = linspace(-2*pi,4*pi,p("N")*3);
  tht1 = linspace(-2*pi,4*pi,p("M")*3)';
  %tht0_n,tht1_n :interpolated grid
  tht0_n = linspace(-2*pi,4*pi,p("num_iter"));
  tht1_n = linspace(-2*pi,4*pi,p("num_iter"));
  % interpolate
  [ri,~] = interpolate_torus(fin_qpos,tht0,tht1,tht0_n,tht1_n,p);

  %% SURF PLOT
  hinter = figure();
  hold on
  grid on
  box on
  axis equal
  % plot orbit -- including periodic orbit
  plot3(L2(1),L2(2),L2(3),'*','MarkerFaceColor','black','MarkerEdgeColor','black','MarkerSize',8);
  view([-24 34])
  xlabel('$x$[-]');
  ylabel('$y$[-]');
  zlabel('$z$[-]');

  %quasi-periodic trajectory
  for i=1:p("N")
    [~,rep] = ode113(@(t,x) fun_cr3bp(t,x,p("mu")),[0 zpo(end)],fin_qpos(p("d")*i-5:p("d")*i),options_ODE);
    hqpt = plot3(rep(:,1),rep(:,2),rep(:,3),"b","Linewidth",0.5);
  end

  %color setting
  CO = zeros(p("num_iter"),p("num_iter"),3);
  CO(:,:,1) = 0.3010.*ones(p("num_iter")); % red
  CO(:,:,2) = 0.5450.*ones(p("num_iter")); % green
  CO(:,:,3) = 0.7.*ones(p("num_iter")); % blue

  %QPT
  hsurf = surf(ri{1,1},ri{2,1},ri{3,1},CO);
  % for visualization
  shading interp
  lightangle(27,36)
  lightangle(27,36)

  legend([hqpt,hsurf],{"Quasi-periodic trajectory","Quasi-periodic tori"},'Location','northeast');
  hold off
end

%% FOCUS ON ONE TORUS
% designate a torus
nm = 5;%1~pacqp("n")
solc_qpos_m = sol_qpos{nm,1};
fin_qpos_m = solc_qpos_m{1,2};

%% DISPLAY CONVERGED INVARIANT CIRCLE
hinvfin = plot_invariant_curve(zpo,fin_qpos_m,p,FR);

%% MANIFOLD SURF PLOT
% calculate the directions of the designated manifolds
del_w_us_inter_full = direction_manifold(fin_qpos_m,p,FR);

% vector arrows
U0 = zeros(p("d")*p("M"),1);

% invariant curve idx
for i=1:p("M") %designate the invariant circle
  k=p("d")*p("N")*(i-1)+1:p("d")*p("N")*i;
  X0 = fin_qpos_m(k);
  del_w_us = del_w_us_inter_full(:,:,i);
  for j=1:p("N")%designate torus point
    % % Create perturbation vector
    % vpert = p("xpert")*norm(X0(p("d")*j-2:p("d")*j))/norm(X0(p("d")*j-5:p("d")*j-3));
    % pert = [ones(3,1).*p("xpert"); ones(3,1).*vpert];

    % Initial points of torus unstable manifolds
    U0(p("d")*i-5:p("d")*i) = X0(p("d")*j-5:p("d")*j)+p("pert").*del_w_us(:,j);
  end
end

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

r1 = sqrt((x+p("mu")).^2+y.^2+z.^2);
r2 = sqrt((x-(1-p("mu"))).^2+y.^2+z.^2);
U = 1/2.*(x.^2+y.^2) + (1-p("mu"))./r1 + p("mu")./r2;
c = 2.*U;
x = reshape(x, [a,b]);
y = reshape(y, [a,b]);
c = reshape(c, [a,b]);

% surf torus
hs = surf(ri{1,1},ri{2,1},ri{3,1},CO);
hold on

% sample trajectories of torus manifolds
% 1th and 2nd points are representives
for i=1:5:p("N")
  [~,X_sample] = ode113(@(t,x) fun_cr3bp(t,x,p("mu")),[0 p("snap_fin_time")+0.1],U0(p("d")*i-5:p("d")*i));
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

