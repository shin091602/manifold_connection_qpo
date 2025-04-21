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
[FR,~,~] = fun_Fourier_interpolation(rho,p);

