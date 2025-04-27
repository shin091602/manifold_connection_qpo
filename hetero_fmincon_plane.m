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
