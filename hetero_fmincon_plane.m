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
p('pert') = 1e-6; % perturbation

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

%% calculate Lyapunov orbits
num_grid = 100;
tspan_1 = linspace(0, 2*lyapunov_init_1(2), num_grid);
[t_corrected_1, x_corrected_1] = ode113(@(t,x) fun_cr3bp(t, x, p('mu')), tspan_1, lyapunov_init_1(3:end), options_ODE);
tspan_2 = linspace(0, 2*lyapunov_init_2(2), num_grid);
[t_corrected_2, x_corrected_2] = ode113(@(t,x) fun_cr3bp(t, x, p('mu')), tspan_2, lyapunov_init_2(3:end), options_ODE);

% plot Lyapunov orbits
figure(1)
plot(x_corrected_1(:, 1), x_corrected_1(:, 2), 'r', 'LineWidth', 2);
hold on
plot(x_corrected_2(:, 1), x_corrected_2(:, 2), 'b', 'LineWidth', 2);
plot(L1(1), L1(2), 'ko', 'MarkerSize', 10, 'MarkerFaceColor', 'k');
plot(L2(1), L2(2), 'ko', 'MarkerSize', 10, 'MarkerFaceColor', 'k');
plot(1-mu, 0, 'ko', 'MarkerSize', 10, 'MarkerFaceColor', 'k');
xlabel('$x$ [-]', 'Interpreter', 'latex', 'FontSize', 24);
ylabel('$y$ [-]', 'Interpreter', 'latex', 'FontSize', 24);
axis equal
grid on
hold off

%% interpolate Lyapunov orbit and monodromy matrix
[~, Y1] = ode113(@(t,x) fun_stm_cr3bp(t,x,p('mu')), tspan_1, [lyapunov_init_1(3:end); reshape(eye(6), 36, 1)], options_ODE);
[~, Y2] = ode113(@(t,x) fun_stm_cr3bp(t,x,p('mu')), tspan_2, [lyapunov_init_2(3:end); reshape(eye(6), 36, 1)], options_ODE);

% Declination angle with respect to L1
theta1 = atan2(Y1(:,1) - L1(1), Y1(:,2) - L1(2));
theta1 = mod(theta1, 2*pi);    % [0,2π] normalization
% sort
[theta1_s, idx1] = sort(theta1);
X1_s = Y1(idx1, :);

% Declination angle with respect to L2
theta2 = atan2(Y2(:,1) - L2(1), Y2(:,2) - L2(2));
theta2 = mod(theta2, 2*pi);    % [0,2π] normalization
% sort
[theta2_s, idx2] = sort(theta2);
X2_s = Y2(idx2, :);

% Generate interpolating functions
interpLyap1 = @(tht_q) cell2mat(arrayfun(@(k) interp1(theta1_s, X1_s(:,k), tht_q, 'pchip'), 1:size(X1_s,2), 'UniformOutput', false) );
interpLyap2 = @(tht_q) cell2mat(arrayfun(@(k) interp1(theta2_s, X2_s(:,k), tht_q, 'pchip'), 1:size(X2_s,2), 'UniformOutput', false) );

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

% Apply perturbation to abitrary point of the Lyapunov orbit
% set the angle of the pointーーーーーーーーーーーーーーーーーーーーーーー
tht_query_1 = pi/6;
tht_query_2 = pi*1.5;
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
U_1 = phi_star_1 * vector_unstable_1;
U_1 = U_1 / norm(U_1);
S_2 = phi_star_2 * vector_stable_2;
S_2 = S_2 / norm(S_2);
U_2 = phi_star_2 * vector_unstable_2;
U_2 = U_2 / norm(U_2);

% create perturbation vector
pert = ones(6, 1) * p('pert');

% Perturb conditions
XS_left_1 = x_star_1' + S_1 .* pert;
XS_right_1 = x_star_1' - S_1 .* pert;
XU_left_1 = x_star_1' + U_1 .* pert;
XU_right_1 = x_star_1' - U_1 .* pert;

XS_left_2 = x_star_2' + S_2 .* pert;
XS_right_2 = x_star_2' - S_2 .* pert;
XU_left_2 = x_star_2' + U_2 .* pert;
XU_right_2 = x_star_2' - U_2 .* pert;

%% calculate the stable and unstable manifolds
tf = 6;
tspan_s = [tf 0];
tspan_u = [0 tf];
[~, xs_left_1] = ode113(@(t,x) fun_cr3bp(t, x, p('mu')), tspan_s, XS_left_1, options_ODE);
[~, xs_right_1] = ode113(@(t,x) fun_cr3bp(t, x, p('mu')), tspan_s, XS_right_1, options_ODE);
[~, xu_left_1] = ode113(@(t,x) fun_cr3bp(t, x, p('mu')), tspan_u, XU_left_1, options_ODE);
[~, xu_right_1] = ode113(@(t,x) fun_cr3bp(t, x, p('mu')), tspan_u, XU_right_1, options_ODE);

[~, xs_left_2] = ode113(@(t,x) fun_cr3bp(t, x, p('mu')), tspan_s, XS_left_2, options_ODE);
[~, xs_right_2] = ode113(@(t,x) fun_cr3bp(t, x, p('mu')), tspan_s, XS_right_2, options_ODE);
[~, xu_left_2] = ode113(@(t,x) fun_cr3bp(t, x, p('mu')), tspan_u, XU_left_2, options_ODE);
[~, xu_right_2] = ode113(@(t,x) fun_cr3bp(t, x, p('mu')), tspan_u, XU_right_2, options_ODE);

%% plot stable and unstable manifolds
figure();
plot(xs_left_1(:, 1), xs_left_1(:, 2), 'r', 'LineWidth', 2);
hold on
plot(xs_right_1(:, 1), xs_right_1(:, 2), 'r', 'LineWidth', 2);
plot(xu_left_1(:, 1), xu_left_1(:, 2), 'b', 'LineWidth', 2);
plot(xu_right_1(:, 1), xu_right_1(:, 2), 'b', 'LineWidth', 2);
plot(xs_left_2(:, 1), xs_left_2(:, 2), 'r', 'LineWidth', 2);
plot(xs_right_2(:, 1), xs_right_2(:, 2), 'r', 'LineWidth', 2);
plot(xu_left_2(:, 1), xu_left_2(:, 2), 'b', 'LineWidth', 2);
plot(xu_right_2(:, 1), xu_right_2(:, 2), 'b', 'LineWidth', 2);
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
axis equal
grid on
hold off

%% calculatate the manifolds to poincare section
% set event function
options_ODE_1   = odeset('RelTol', 1e-13, 'AbsTol', 1e-13, 'Events', @(t,x) odestop_hetero_1(t,x,mu));
options_ODE_2   = odeset('RelTol', 1e-13, 'AbsTol', 1e-13, 'Events', @(t,x) odestop_hetero_2(t,x,mu));
% calculate the manifolds
[~, xs_left_1, ~, xes_left_1, ~] = ode113(@(t,x) fun_cr3bp(t, x, p('mu')), tspan_s, XS_left_1, options_ODE_1);
[~, xs_right_1, ~, xes_right_1, ~] = ode113(@(t,x) fun_cr3bp(t, x, p('mu')), tspan_s, XS_right_1, options_ODE_1);
[~, xu_left_1, ~, xeu_left_1, ~] = ode113(@(t,x) fun_cr3bp(t, x, p('mu')), tspan_u, XU_left_1, options_ODE_1);
[~, xu_right_1, ~, xeu_right_1, ~] = ode113(@(t,x) fun_cr3bp(t, x, p('mu')), tspan_u, XU_right_1, options_ODE_1);

[~, xs_left_2, ~, xes_left_2, ~] = ode113(@(t,x) fun_cr3bp(t, x, p('mu')), tspan_s, XS_left_2, options_ODE_2);
[~, xs_right_2, ~, xes_right_2, ~] = ode113(@(t,x) fun_cr3bp(t, x, p('mu')), tspan_s, XS_right_2, options_ODE_2);
[~, xu_left_2, ~ , xeu_left_2 ,~] = ode113(@(t,x) fun_cr3bp(t,x,p('mu')), tspan_u,XU_left_2 ,options_ODE_2);
[~, xu_right_2 ,~ , xeu_right_2 ,~] = ode113(@(t,x) fun_cr3bp(t,x,p('mu')), tspan_u,XU_right_2 ,options_ODE_2);