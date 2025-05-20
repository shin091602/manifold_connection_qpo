% CR3BP
% For plot heterocrinic connections in the CR3BP planetary system

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code uses fmincon to compute heteroclinic connections
% that exist between Lyapunov orbits in a plane.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

%% initialization
[mu ,~,~,~] = parameter(2); %mu = 1.21536E-02; % Earth-Moon
[L1, L2, ~, ~, ~] = librationPoints(mu); %Lagrange points

% ODE options
options_ODE = odeset('RelTol',1e-13,'AbsTol',1e-13);
options_ODE_1   = odeset('RelTol', 1e-13, 'AbsTol', 1e-13, 'Events', @(t,x) odestop_hetero_1(t,x,mu));
options_ODE_2   = odeset('RelTol', 1e-13, 'AbsTol', 1e-13, 'Events', @(t,x) odestop_hetero_2(t,x,mu));

%Dictionary of the variables
% constant values for the calculation are defined here

p = dictionary();

p("Me") = 5.9724e+24;%mass earth
p("M1") = 7.3458e+22;%mass moon
p("mu") = 1.21536E-02;%mu--EM
p("G") = 6.67300000000000e-11;%gravitational constant
p("chara_length_CR3BP") = 3.844e+8;%[m]
p("chara_mass_CR3BP") = 5.9724e+24;%[kg]
p("chara_time_CR3BP") = 3.7748e+5;%[s]
p("N_CR3BP_SE") = 1.6645e-05;

p("d") = 6; % dimension of variables

% target Jacobi constant
p("C_xn") = 3.125;
p("C_error") = 1e-13;

% setting for continuation
p('iteration_max') = 1000;
p('threshold') = 1e-10;
p('count_max') = 8000;

% setting for propagation of the manifolds
p('pert') = 1e-8; % perturbation

%% natural parameter continuation
delta = 1e-4; % step size
count = 0; % count of the iteration

lyapunov_init_1 = zeros(8, 1);
lyapunov_init_2 = zeros(8, 1);

x0_1 = [0.836900082907655, 0, 0, 0, 1.770874936727959e-06, 0];
t0_1 = 1.3458;

x0_2 = [1.155693450711950, 0, 0, 0, 1.666612527554842e-06, 0];
t0_2 = 1.68665;

%% Loop for L1-Lyapunov orbit
while 1
    count = count + 1;

    % differential correction
    for iteration = 1:p('iteration_max')
        [x_n_1, t_n_1, C] = fun_differential_correction_cr3bp(x0_1, t0_1, p('mu'),options_ODE);

        tspan = [0 2*t_n_1];
        [t_corrected_1, x_corrected_1] = ode113(@(t,x) fun_cr3bp(t, x, p('mu')), tspan, x_n_1, options_ODE);

        x_error_1 = norm(x_corrected_1(end, :) - x_corrected_1(1, :));

        if x_error_1 < p('threshold')
            break;
        end

        if x_error_1 > 1e+3
            disp('calculation diverged')
            return;
        end

        if iteration == p('iteration_max')
            disp('reach iteration max')
        end

        x0_1 = x_n_1;
        t0_1 = t_n_1;
    end

    if C < p("C_xn")
        x0_1(1) = x0_1(1) - delta;
        delta = delta / 5;
        disp(strcat('delta changed : count = ', num2str(count), ' delta = ', num2str(delta)));
    end

    if (count >= 100)&&(mod(count,100)==0)
        disp(strcat('count = ', num2str(count), ' delta = ', num2str(delta)));
    end

    if abs(C - p("C_xn")) < p("C_error")
        lyapunov_init_1 = [C; t_n_1; x_n_1(:,1)];
        disp('L1 Lyapunov orbit found');
        break
    end

    x0_1(1) = x0_1(1) + delta;

    if count == p('count_max') - 1
        break
    end
end

%% Loop for L2-Lyapunov orbit
count = 0;
delta = 1e-4; % step size

while 1
    count = count + 1;

    % differential correction
    for iteration = 1:p('iteration_max')
        [x_n_2, t_n_2, C] = fun_differential_correction_cr3bp(x0_2, t0_2, p('mu'),options_ODE);

        tspan = [0 2*t_n_2];
        [t_corrected_2, x_corrected_2] = ode113(@(t,x) fun_cr3bp(t, x, p('mu')), tspan, x_n_2, options_ODE);

        x_error_2 = norm(x_corrected_2(end, :) - x_corrected_2(1, :));

        if x_error_2 < p('threshold')
            break;
        end

        if x_error_2 > 1e+3
            disp('calculation diverged')
            return;
        end

        if iteration == p('iteration_max')
            disp('reach iteration max')
        end

        x0_2 = x_n_2;
        t0_2 = t_n_2;
    end

    if C < p("C_xn")
        x0_2(1) = x0_2(1) - delta;
        delta = delta / 5;
        disp(strcat('delta changed : count = ', num2str(count), ' delta = ', num2str(delta)));
    end

    if (count >= 100)&&(mod(count,100)==0)
        disp(strcat('count = ', num2str(count), ' delta = ', num2str(delta)));
    end

    if abs(C - p("C_xn")) < p("C_error")
        lyapunov_init_2 = [C; t_n_2; x_n_2(:,1)];
        disp('L2 Lyapunov orbit found');
        break
    end

    x0_2(1) = x0_2(1) + delta;

    if count == p('count_max') - 1
        break
    end
end

%% check the mainifolds before interpolation
N = 100;
[~, XS_right_1,~,~] = fun_manifold_cr3bp(mu,lyapunov_init_1(3:end),lyapunov_init_1(2),N,p('pert'),options_ODE_1);
[~,~,XU_right_2,~] = fun_manifold_cr3bp(mu,lyapunov_init_2(3:end),lyapunov_init_2(2),N,p('pert'),options_ODE_2);
tf = 15;
tspan_s = [tf 0];
tspan_u = [0 tf];
color      = jet;
N_lim = linspace(0, N, size(color,1));
Interp_N   = griddedInterpolant(N_lim, color);
figure();
hold on
for i = 1:N
    [~,ys_1,~,~,~] = ode113(@(t,x) fun_cr3bp(t,x,p('mu')), tspan_s, XS_right_1(:,i), options_ODE_1);
    plot(ys_1(:,1), ys_1(:,2), 'Color', Interp_N(i),'HandleVisibility','off');
    [~,yu_2,~,~,~] = ode113(@(t,x) fun_cr3bp(t,x,p('mu')), tspan_u, XU_right_2(:,i), options_ODE_2);
    plot(yu_2(:,1), yu_2(:,2), 'Color', Interp_N(i),'HandleVisibility','off');
end
plot(x_corrected_1(:, 1), x_corrected_1(:, 2), 'k', 'LineWidth', 2,'HandleVisibility','off');
plot(x_corrected_2(:, 1), x_corrected_2(:, 2), 'k', 'LineWidth', 2,'HandleVisibility','off');
plot(L1(1), L1(2), 'k*', 'MarkerSize', 10, 'MarkerFaceColor', 'k','DisplayName','$L_1$');
plot(L2(1), L2(2), 'k+', 'MarkerSize', 10, 'MarkerFaceColor', 'k','DisplayName','$L_2$');
plot(1-mu, 0, 'ko', 'MarkerSize', 10, 'MarkerFaceColor', 'k','DisplayName','Moon');
xlabel('$x$ [-]', 'Interpreter', 'latex', 'FontSize', 24);
ylabel('$y$ [-]', 'Interpreter', 'latex', 'FontSize', 24);
legend('Location','best');
xlim([0.7 1.3]);
ylim([-0.2, 0.2]);
grid on
title('Manifolds before interpolation');
hold off

%% check the poincare section before interpolation
figure();
hold on
N=500;
color      = jet;
N_lim = linspace(0, N, size(color,1));
Interp_N   = griddedInterpolant(N_lim, color);
[~, XS_right_1,~,~] = fun_manifold_cr3bp(mu,lyapunov_init_1(3:end),lyapunov_init_1(2),N,p('pert'),options_ODE_1);
[~,~,XU_right_2,~] = fun_manifold_cr3bp(mu,lyapunov_init_2(3:end),lyapunov_init_2(2),N,p('pert'),options_ODE_2);
for i = 1:N
    [~,~,~,ye1,~] = ode113(@(t,x) fun_cr3bp(t,x,p('mu')), tspan_s, XS_right_1(:,i), options_ODE_1);
    if ye1
        if ye1(2) > 0.1
            plot3(ye1(2), ye1(5), i, '+', 'Color', Interp_N(i), 'MarkerSize', 5);
        end
    end
    [~,~,~,ye2,~] = ode113(@(t,x) fun_cr3bp(t,x,p('mu')), tspan_u, XU_right_2(:,i), options_ODE_2);
    if ye2
        if ye2(2) > 0.1
            plot3(ye2(2), ye2(5), i, '*', 'Color', Interp_N(i), 'MarkerSize', 5);
        end
    end
end
hL1 = plot(nan, nan, '+', 'MarkerSize', 5, 'Color', [0 0 0]);  % L1 マーカー
hL2 = plot(nan, nan, '*', 'MarkerSize', 5, 'Color', [0 0 0]);  % L2 マーカー
% 凡例を出力
legend([hL1, hL2], {'L1 stable manifold', 'L2 unstable manifold'}, 'Location', 'best');
xlabel('$y$ [-]', 'Interpreter', 'latex', 'FontSize', 24);
ylabel('$dy$ [-]', 'Interpreter', 'latex', 'FontSize', 24);
zlabel('$number$ [-]', 'Interpreter', 'latex', 'FontSize', 24);
title('Poincare section before interpolation');
grid on
hold off
%% calculate Lyapunov orbits
num_grid = 1000;
tspan_1 = linspace(0, 2*lyapunov_init_1(2), num_grid);
[t_corrected_1, x_corrected_1] = ode113(@(t,x) fun_cr3bp(t, x, p('mu')), tspan_1, lyapunov_init_1(3:end), options_ODE);
tspan_2 = linspace(0, 2*lyapunov_init_2(2), num_grid);
[t_corrected_2, x_corrected_2] = ode113(@(t,x) fun_cr3bp(t, x, p('mu')), tspan_2, lyapunov_init_2(3:end), options_ODE);

% plot Lyapunov orbits
figure();
hold on
plot(x_corrected_1(:, 1), x_corrected_1(:, 2), 'r', 'LineWidth', 2);
plot(x_corrected_1(1,1), x_corrected_1(1,2), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
plot(x_corrected_2(:, 1), x_corrected_2(:, 2), 'b', 'LineWidth', 2);
plot(x_corrected_2(1,1), x_corrected_2(1,2), 'bo', 'MarkerSize', 10, 'MarkerFaceColor', 'b');
plot(L1(1), L1(2), 'ko', 'MarkerSize', 10, 'MarkerFaceColor', 'k');
plot(L2(1), L2(2), 'ko', 'MarkerSize', 10, 'MarkerFaceColor', 'k');
plot(1-mu, 0, 'ko', 'MarkerSize', 10, 'MarkerFaceColor', 'k');
xlabel('$x$ [-]', 'Interpreter', 'latex', 'FontSize', 24);
ylabel('$y$ [-]', 'Interpreter', 'latex', 'FontSize', 24);
axis equal
grid on
hold off

%% interpolate Lyapunov orbit and monodromy matrix
[t1, Y1] = ode113(@(t,x) fun_stm_cr3bp(t,x,p('mu')), tspan_1, [lyapunov_init_1(3:end); reshape(eye(6), 36, 1)], options_ODE);
[t2, Y2] = ode113(@(t,x) fun_stm_cr3bp(t,x,p('mu')), tspan_2, [lyapunov_init_2(3:end); reshape(eye(6), 36, 1)], options_ODE);

% Declination angle with respect to L1
theta1 = atan2(Y1(:,2) - L1(2), Y1(:,1) - L1(1));
figure();
plot(theta1, 'r');
grid on
title('theta1');
theta1 = mod(theta1, 2*pi);    % [0,2π] normalization
figure();
plot(theta1, 'r');
grid on
title('theta1_mod');
% sort
[theta1_s, idx1] = sort(theta1);
figure();
plot(theta1_s, 'r');
grid on
title('theta1_sort');
X1_s = Y1(idx1, :);

% Declination angle with respect to L2
theta2 = atan2(Y2(:,2) - L2(2), Y2(:,1) - L2(1));
theta2 = mod(theta2, 2*pi);    % [0,2π] normalization
% sort
[theta2_s, idx2] = sort(theta2);
X2_s = Y2(idx2, :);

% Generate interpolating functions
interpLyap1 = @(tht_q) cell2mat(arrayfun(@(k) interp1(theta1_s, X1_s(:,k), tht_q, 'pchip'), 1:size(X1_s,2), 'UniformOutput', false) );
interpLyap2 = @(tht_q) cell2mat(arrayfun(@(k) interp1(theta2_s, X2_s(:,k), tht_q, 'pchip'), 1:size(X2_s,2), 'UniformOutput', false) );
%% check theta
theta_range = linspace(0, 2*pi, 1000);
color      = jet;
theta_lim = linspace(0, 2*pi, size(color,1));
Interp_theta   = griddedInterpolant(theta_lim, color);
figure();
hold on
for i = 1:length(theta_range)
    % L1
    X1 = interpLyap1(theta_range(i));
    x1 = X1(1:6);
    plot(x1(1), x1(2), '+','color', Interp_theta(theta_range(i)), 'MarkerSize', 5,'HandleVisibility','off');
    if i == 1
        plot(x1(1), x1(2), 'o', 'color', Interp_theta(theta_range(i)), 'MarkerSize', 20,'DisplayName', 'theta = 0');
    end

    % L2
    X2 = interpLyap2(theta_range(i));
    x2 = X2(1:6);
    plot(x2(1), x2(2),'+', 'color', Interp_theta(theta_range(i)), 'MarkerSize', 5,'HandleVisibility','off');
    if i == 1
        plot(x2(1), x2(2), 'o', 'color', Interp_theta(theta_range(i)), 'MarkerSize', 20,'DisplayName', 'theta = 0');
    end
end
grid on
legend('Location','best');
hold off
%% use the interpolating functions
% monodromy matrix
M_1 = reshape(Y1(end, 7:end), 6, 6);
M_2 = reshape(Y2(end, 7:end), 6, 6);
% Eigenvectors and eigenvalues analysis
[V1, D1] = eig(M_1);
D1 = diag(D1);
[~, idx_s1] = min(abs(D1)); % stable
[~, idx_u1] = max(abs(D1)); % unstable
[V2, D2] = eig(M_2);
D2 = diag(D2);
[~, idx_s2] = min(abs(D2)); % stable
[~, idx_u2] = max(abs(D2)); % unstable

% Stable and unstable eigenvectors
vector_stable_1 = V1(:, idx_s1);
if vector_stable_1(1) < 0
    vector_stable_1 = -vector_stable_1;
end
vector_unstable_1 = V1(:, idx_u1);
if vector_unstable_1(1) < 0
    vector_unstable_1 = -vector_unstable_1;
end

vector_stable_2 = V2(:, idx_s2);
if vector_stable_2(1) < 0
    vector_stable_2 = -vector_stable_2;
end
vector_unstable_2 = V2(:, idx_u2);
if vector_unstable_2(1) < 0
    vector_unstable_2 = -vector_unstable_2;
end

% %% Apply perturbation to abitrary point of the Lyapunov orbit
% % set the angle of the pointーーーーーーーーーーーーーーーーーーーーーーー
% tht_query_1 = 0.52501570;
% tht_query_2 = 4.64640599;
% X1_resampled = interpLyap1(tht_query_1);
% X2_resampled = interpLyap2(tht_query_2);
% % Grab state at the fixed point
% x_star_1 = X1_resampled(1:6);
% x_star_2 = X2_resampled(1:6);
% % Grab state transition matrix at the fixed point
% phi_star_1 = reshape(X1_resampled(7:end), 6, 6);
% phi_star_2 = reshape(X2_resampled(7:end), 6, 6);
%
% % Map stable and unstable vectors forward
% S_1 = phi_star_1 * vector_stable_1;
% S_1 = S_1 / norm(S_1);
% U_1 = phi_star_1 * vector_unstable_1;
% U_1 = U_1 / norm(U_1);
% S_2 = phi_star_2 * vector_stable_2;
% S_2 = S_2 / norm(S_2);
% U_2 = phi_star_2 * vector_unstable_2;
% U_2 = U_2 / norm(U_2);
%
% % create perturbation vector
% pert = ones(6, 1) * p('pert');
%
% % Perturb conditions
% XS_left_1 = x_star_1' + S_1 .* pert;
% XS_right_1 = x_star_1' - S_1 .* pert;
% XU_left_1 = x_star_1' + U_1 .* pert;
% XU_right_1 = x_star_1' - U_1 .* pert;
%
% XS_left_2 = x_star_2' + S_2 .* pert;
% XS_right_2 = x_star_2' - S_2 .* pert;
% XU_left_2 = x_star_2' + U_2 .* pert;
% XU_right_2 = x_star_2' - U_2 .* pert;
%
% %% calculate the stable and unstable manifolds
tf = 15;
tspan_s = [tf 0];
tspan_u = [0 tf];
%% search the initial point of the heteroclinic connection
theta_range = linspace(0, 2*pi, 5000);
color      = jet;
theta_lim = linspace(0, 2*pi, size(color,1));
Interp_theta   = griddedInterpolant(theta_lim, color);
figure();
hold on
for i = 1:length(theta_range)
    % L1
    X1 = interpLyap1(theta_range(i));
    x1 = X1(1:6);
    phi1 = reshape(X1(7:end), 6, 6);
    S = phi1 * vector_stable_1;
    S = S / norm(S);
    XS = x1' + p('pert') * S;
    [~, ~, ~, YE1, ~] = ode113(@(t,x) fun_cr3bp(t,x,p('mu')), tspan_s, XS, options_ODE_1);
    if YE1
        if YE1(2) > 0.1
            plot3(YE1(2), YE1(5), theta_range(i), '+', 'Color', Interp_theta(theta_range(i)), 'MarkerSize', 5);
        end
    end

    % L2
    X2 = interpLyap2(theta_range(i));
    x2 = X2(1:6);
    phi2 = reshape(X2(7:end), 6, 6);
    U = phi2 * vector_stable_2;
    U = U / norm(U);
    XU = x2' - p('pert') * U;
    [~, ~, ~, YE2, ~] = ode113(@(t,x) fun_cr3bp(t,x,p('mu')), tspan_u, XU, options_ODE_2);
    if YE2
        if YE2(2) > 0.1
            plot3(YE2(2), YE2(5), theta_range(i), '*', 'Color', Interp_theta(theta_range(i)), 'MarkerSize', 5);
        end
    end
end
hL1 = plot(nan, nan, '+', 'MarkerSize', 5, 'Color', [0 0 0]);  % L1 マーカー
hL2 = plot(nan, nan, '*', 'MarkerSize', 5, 'Color', [0 0 0]);  % L2 マーカー
% 凡例を出力
legend([hL1, hL2], {'L1 stable manifold', 'L2 unstable manifold'}, 'Location', 'best');
colormap jet;
c = colorbar;
ylabel(c, 'theta [-]', 'FontSize', 15);
caxis([0 2*pi]);
c.Ticks = linspace(0, 2*pi, 20);
xlim([-0.2 0.2]);
ylim([-0.5 0.5]);
grid on
hold off

%% check where manifolds come from
theta_range = linspace(0, 2*pi, 250);
color      = jet;
theta_lim = linspace(0, 2*pi, size(color,1));
Interp_theta   = griddedInterpolant(theta_lim, color);
figure();
hold on
for i = 1:length(theta_range)
    % L1
    X1 = interpLyap1(theta_range(i));
    x1 = X1(1:6);
    phi1 = reshape(X1(7:end), 6, 6);
    S = phi1 * vector_stable_1;
    S = S / norm(S);
    XS = x1' + p('pert') * S;
    [~, ys_1, ~, YE1, ~] = ode113(@(t,x) fun_cr3bp(t,x,p('mu')), tspan_s, XS, options_ODE_1);
    if YE1
        plot(ys_1(:,1), ys_1(:,2), 'Color', Interp_theta(theta_range(i)));
    end

    % L2
    X2 = interpLyap2(theta_range(i));
    x2 = X2(1:6);
    phi2 = reshape(X2(7:end), 6, 6);
    U = phi2 * vector_stable_2;
    U = U / norm(U);
    XU = x2' - p('pert') * U;
    [~, yu_2, ~, YE2, ~] = ode113(@(t,x) fun_cr3bp(t,x,p('mu')), tspan_u, XU, options_ODE_2);
    if YE2
        plot(yu_2(:,1), yu_2(:,2),'Color', Interp_theta(theta_range(i)));
    end
end
colormap jet;
c = colorbar;
ylabel(c, 'theta [-]', 'FontSize', 15);
caxis([0 2*pi]);
c.Ticks = linspace(0, 2*pi, 20);
grid on
hold off

%% use fsolve to find the intersection
% 初期値（お好みで調整）
th0 = [5.1507; 5.0212];

optsFS = optimoptions('fsolve', ...
    'Display','iter', ...
    'TolFun',1e-12, ...
    'TolX',1e-12, ...
    'Algorithm','levenberg-marquardt');  % m>=n の系に適したアルゴリズム

% handle function
systemFun = @(th) hetero_system(th, interpLyap1, interpLyap2,vector_unstable_1, vector_stable_2,p, tspan_u, tspan_s, options_ODE_1, options_ODE_2);

% solve
[th_sol, Fval, ~] = fsolve(systemFun, th0, optsFS);
fprintf('θ₁* = %.8f, θ₂* = %.8f, residual norm = %.3e\n', th_sol(1), th_sol(2), norm(Fval));
% Apply perturbation to abitrary point of the Lyapunov orbit
% set the angle of the pointーーーーーーーーーーーーーーーーーーーーーーー
tht_query_1 = th_sol(1);
tht_query_2 = th_sol(2);
X1_resampled = interpLyap1(tht_query_1);
X2_resampled = interpLyap2(tht_query_2);
% Grab state at the fixed point
x_star_1 = X1_resampled(1:6);
x_star_2 = X2_resampled(1:6);
% Grab state transition matrix at the fixed point
phi_star_1 = reshape(X1_resampled(7:end), 6, 6);
phi_star_2 = reshape(X2_resampled(7:end), 6, 6);

% Map stable and unstable vectors forward
S_1 = phi_star_1 * vector_stable_1;
S_1 = S_1 / norm(S_1);
U_2 = phi_star_2 * vector_unstable_2;
U_2 = U_2 / norm(U_2);

% create perturbation vector
pert = ones(6, 1) * p('pert');

% Perturb conditions
XS_right_1 = x_star_1' + S_1 .* pert;
XU_left_2 = x_star_2' - U_2 .* pert;
% calculate the stable and unstable manifolds
tf = 30;
tspan_s = [tf 0];
tspan_u = [0 tf];

% calculatate the manifolds to poincare section
% calculate the manifolds
[~, xs_right_1, ~, xes_right_1, ~] = ode113(@(t,x) fun_cr3bp(t, x, p('mu')), tspan_s, XS_right_1, options_ODE_1);
[~, xu_left_2, ~ , xeu_left_2 ,~] = ode113(@(t,x) fun_cr3bp(t,x,p('mu')), tspan_u,XU_left_2 ,options_ODE_2);

% plot stable and unstable manifolds
figure();
hold on
plot(xs_right_1(:, 1), xs_right_1(:, 2), 'r', 'LineWidth', 2);
plot(xu_left_2(:, 1), xu_left_2(:, 2), 'b', 'LineWidth', 2);
plot(x_corrected_1(:, 1), x_corrected_1(:, 2), 'k', 'LineWidth', 2);
plot(x_corrected_2(:, 1), x_corrected_2(:, 2), 'k', 'LineWidth', 2);
plot(L1(1), L1(2), 'ko', 'MarkerSize', 10, 'MarkerFaceColor', 'k');
plot(L2(1), L2(2), 'ko', 'MarkerSize', 10, 'MarkerFaceColor', 'k');
plot(x_star_1(1), x_star_1(2), 'ko', 'MarkerSize', 10, 'MarkerFaceColor', 'm');
plot(x_star_2(1), x_star_2(2), 'ko', 'MarkerSize', 10, 'MarkerFaceColor', 'm');
plot(1-mu, 0, 'ko', 'MarkerSize', 10, 'MarkerFaceColor', 'k');
xlabel('$x$ [-]', 'Interpreter', 'latex', 'FontSize', 24);
ylabel('$y$ [-]', 'Interpreter', 'latex', 'FontSize', 24);
xlim([0.7 1.3]);
ylim([-0.2, 0.2]);
title('fsolve')
axis equal
grid on
hold off

%% fmincon
th0 = [5.1507; 5.0212];

% handle function
systemFun = @(th) hetero_system(th, interpLyap1, interpLyap2,vector_unstable_1, vector_stable_2,p, tspan_u, tspan_s, options_ODE_1, options_ODE_2);


% fmincon 用のオプション設定
optsFC = optimoptions('fmincon', ...
    'Display', 'iter', ...           % 反復過程を表示
    'Algorithm', 'interior-point', ... % 制約付き最適化に適した手法
    'TolFun', 1e-12, ...
    'TolX',   1e-12, ...
    'MaxIterations', 400);

% 目的関数の定義：残差ノルムの二乗を最小化
objFun = @(th) sum( systemFun(th).^2 );

% 制約がなければ下記のように空配列で与えられる
A     = [];  b     = [];
Aeq   = [];  beq   = [];
lb    = [];  ub    = [];
nonlcon = [];

% 最適化実行
[th_sol, fval_sq, exitflag, output] = fmincon( ...
    objFun, th0, A, b, Aeq, beq, lb, ub, nonlcon, optsFC);

% fval_sq は「残差ノルムの二乗」の最小値なので、実際のノルムは
resnorm = sqrt(fval_sq);
fprintf('θ₁* = %.8f, θ₂* = %.8f, residual norm = %.3e\n', ...
    th_sol(1), th_sol(2), resnorm);

% set the angle of the pointーーーーーーーーーーーーーーーーーーーーーーー
tht_query_1 = th_sol(1);
tht_query_2 = th_sol(2);
X1_resampled = interpLyap1(tht_query_1);
X2_resampled = interpLyap2(tht_query_2);
% Grab state at the fixed point
x_star_1 = X1_resampled(1:6);
x_star_2 = X2_resampled(1:6);
% Grab state transition matrix at the fixed point
phi_star_1 = reshape(X1_resampled(7:end), 6, 6);
phi_star_2 = reshape(X2_resampled(7:end), 6, 6);

% Map stable and unstable vectors forward
S_1 = phi_star_1 * vector_stable_1;
S_1 = S_1 / norm(S_1);
U_2 = phi_star_2 * vector_unstable_2;
U_2 = U_2 / norm(U_2);

% create perturbation vector
pert = ones(6, 1) * p('pert');

% Perturb conditions
XS_right_1 = x_star_1' + S_1 .* pert;
XU_left_2 = x_star_2' - U_2 .* pert;
% calculate the stable and unstable manifolds
tf = 30;
tspan_s = [tf 0];
tspan_u = [0 tf];

% calculatate the manifolds to poincare section
% calculate the manifolds
[~, xs_right_1, ~, xes_right_1, ~] = ode113(@(t,x) fun_cr3bp(t, x, p('mu')), tspan_s, XS_right_1, options_ODE_1);
[~, xu_left_2, ~ , xeu_left_2 ,~] = ode113(@(t,x) fun_cr3bp(t,x,p('mu')), tspan_u,XU_left_2 ,options_ODE_2);

% plot stable and unstable manifolds
figure();
hold on
plot(xs_right_1(:, 1), xs_right_1(:, 2), 'r', 'LineWidth', 2);
plot(xu_left_2(:, 1), xu_left_2(:, 2), 'b', 'LineWidth', 2);
plot(x_corrected_1(:, 1), x_corrected_1(:, 2), 'k', 'LineWidth', 2);
plot(x_corrected_2(:, 1), x_corrected_2(:, 2), 'k', 'LineWidth', 2);
plot(L1(1), L1(2), 'ko', 'MarkerSize', 10, 'MarkerFaceColor', 'k');
plot(L2(1), L2(2), 'ko', 'MarkerSize', 10, 'MarkerFaceColor', 'k');
plot(x_star_1(1), x_star_1(2), 'ko', 'MarkerSize', 10, 'MarkerFaceColor', 'm');
plot(x_star_2(1), x_star_2(2), 'ko', 'MarkerSize', 10, 'MarkerFaceColor', 'm');
plot(1-mu, 0, 'ko', 'MarkerSize', 10, 'MarkerFaceColor', 'k');
xlabel('$x$ [-]', 'Interpreter', 'latex', 'FontSize', 24);
ylabel('$y$ [-]', 'Interpreter', 'latex', 'FontSize', 24);
xlim([0.7 1.3]);
ylim([-0.2, 0.2]);
title('fmincon')
axis equal
grid on
hold off


%% fsolve version2
% 初期値（お好みで調整）
th0 = [5.1507; 5.0212];

% fsolve オプション
optsFS = optimoptions('fsolve', ...
    'Display','iter', ...
    'TolFun',1e-12, ...
    'TolX',1e-12);

% システム関数のハンドル
systemFun = @(th) hetero_system(th, interpLyap1, interpLyap2, ...
    vector_unstable_1, vector_stable_2, ...
    p, tspan_u, tspan_s, ...
    options_ODE_1, options_ODE_2);

% 解を求める
[th_sol, Fval, exitflag] = fsolve(systemFun, th0, optsFS);
fprintf('θ₁* = %.8f, θ₂* = %.8f, residual norm = %.3e\n', ...
    th_sol(1), th_sol(2), norm(Fval));

tht_query_1 = th_sol(1);
tht_query_2 = th_sol(2);
X1_resampled = interpLyap1(tht_query_1);
X2_resampled = interpLyap2(tht_query_2);
% Grab state at the fixed point
x_star_1 = X1_resampled(1:6);
x_star_2 = X2_resampled(1:6);
% Grab state transition matrix at the fixed point
phi_star_1 = reshape(X1_resampled(7:end), 6, 6);
phi_star_2 = reshape(X2_resampled(7:end), 6, 6);

% Map stable and unstable vectors forward
S_1 = phi_star_1 * vector_stable_1;
S_1 = S_1 / norm(S_1);
U_2 = phi_star_2 * vector_unstable_2;
U_2 = U_2 / norm(U_2);

% create perturbation vector
pert = ones(6, 1) * p('pert');

% Perturb conditions
XS_right_1 = x_star_1' + S_1 .* pert;
XU_left_2 = x_star_2' - U_2 .* pert;
% calculate the stable and unstable manifolds
tf = 30;
tspan_s = [tf 0];
tspan_u = [0 tf];

% calculatate the manifolds to poincare section
% calculate the manifolds
[~, xs_right_1, ~, xes_right_1, ~] = ode113(@(t,x) fun_cr3bp(t, x, p('mu')), tspan_s, XS_right_1, options_ODE_1);
[~, xu_left_2, ~ , xeu_left_2 ,~] = ode113(@(t,x) fun_cr3bp(t,x,p('mu')), tspan_u,XU_left_2 ,options_ODE_2);

% plot stable and unstable manifolds
figure();
hold on
plot(xs_right_1(:, 1), xs_right_1(:, 2), 'r', 'LineWidth', 2);
plot(xu_left_2(:, 1), xu_left_2(:, 2), 'b', 'LineWidth', 2);
plot(x_corrected_1(:, 1), x_corrected_1(:, 2), 'k', 'LineWidth', 2);
plot(x_corrected_2(:, 1), x_corrected_2(:, 2), 'k', 'LineWidth', 2);
plot(L1(1), L1(2), 'ko', 'MarkerSize', 10, 'MarkerFaceColor', 'k');
plot(L2(1), L2(2), 'ko', 'MarkerSize', 10, 'MarkerFaceColor', 'k');
plot(x_star_1(1), x_star_1(2), 'ko', 'MarkerSize', 10, 'MarkerFaceColor', 'm');
plot(x_star_2(1), x_star_2(2), 'ko', 'MarkerSize', 10, 'MarkerFaceColor', 'm');
plot(1-mu, 0, 'ko', 'MarkerSize', 10, 'MarkerFaceColor', 'k');
xlabel('$x$ [-]', 'Interpreter', 'latex', 'FontSize', 24);
ylabel('$y$ [-]', 'Interpreter', 'latex', 'FontSize', 24);
xlim([0.7 1.3]);
ylim([-0.2, 0.2]);
title('fsolve_version2')
axis equal
grid on
hold off

%%
% 福岡空港と羽田空港の緯度・経度（度単位）
lat1 = 33 + 35/60 + 4/3600;    % 福岡空港 ARP 緯度 [°]
lon1 = 130 + 27/60 + 6/3600;   % 福岡空港 ARP 経度 [°]
lat2 = 35 + 33/60 + 12/3600;   % 羽田空港 ARP 緯度 [°]
lon2 = 139 + 46/60 + 52/3600;  % 羽田空港 ARP 経度 [°]

% 度→ラジアン変換（ASCII変数名）
phi1   = deg2rad(lat1);
lam1   = deg2rad(lon1);
phi2   = deg2rad(lat2);
lam2   = deg2rad(lon2);

% 経度差
deltaLam = lam2 - lam1;

% 方位角計算（北を0°、時計回り）
thetaRad = atan2( sin(deltaLam) * cos(phi2), ...
    cos(phi1) * sin(phi2) - sin(phi1) * cos(phi2) * cos(deltaLam) );

% 0–360° に正規化
bearingDeg = mod(rad2deg(thetaRad) + 360, 360);

fprintf('福岡ARP→羽田ARP の方位角: %.3f°\n', bearingDeg);
