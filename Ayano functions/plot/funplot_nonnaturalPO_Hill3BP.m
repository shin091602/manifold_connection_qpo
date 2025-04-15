function [f1, t_ini, x_ini, initial_new, period_new] = funplot_nonnaturalPO_Hill3BP(x0, t0)
% t0 : T
% x0 : initial

%% differential correction
iteration_max = 100;
options_ODE = odeset('RelTol',1e-14, 'AbsTol',1e-14);

if (abs(x0(3)) < 1e-7) && (abs(x0(6)) < 1e-7) % Lyapunov
    threshold = 1e-12;
    for s = 1 : iteration_max
        [initial_new, period_new, ~, ~] = fun_2n_differential_correction_Hill3BP_NPC_x0fixed(x0, t0, options_ODE, Q, xt);

        tspan = [0 period_new];
        [t_ini, x_ini] = ode113(@(t, xp) fun_ode_2n_Hill3BP(t, xp, Q, xt), tspan, initial_new, options_ODE);

        error = norm(x_ini(end,:) - x_ini(1,:));
        if error < threshold % 完全に0だと数値誤差の関係で一生収束しないので
            break;
        end

        if error > 1e+3
            disp('calculation diverged');
            return;
        end

        if s == iteration_max
            disp('do not finish');
            return;
        end
    end
elseif (abs(x0(3)) > 1e-4) && (abs(x0(6)) < 1e-7) % halo
    threshold = 1e-12;
    for s = 1 : iteration_max
        [initial_new, period_new, ~, ~] = fun_2n_differential_correction_Hill3BP_NPC_x0fixed(x0, t0, options_ODE, Q, xt);

        tspan = [0 period_new];
        [t_ini, x_ini] = ode113(@(t, xp) fun_ode_2n_Hill3BP(t, xp, Q, xt), tspan, initial_new, options_ODE);

        error = norm(x_ini(end,:) - x_ini(1,:));
        if error < threshold % 完全に0だと数値誤差の関係で一生収束しないので
            break;
        end

        if error > 1e+3
            disp('calculation diverged');
            return;
        end

        if s == iteration_max
            disp('do not finish');
            return;
        end
    end
end
x_ini = x_ini';

%% plot
f1=figure;
axis equal
hold on
plot3(x_ini(1,:), x_ini(2,:), x_ini(3,:));
hold off;
end