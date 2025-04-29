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
        [x_n_1, t_n_1, C] = fun_differential_correction_cr3bp(x0_1, t0_1, mu,options_ODE);

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


%%
