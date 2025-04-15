function [f1, t_ini, x_ini, initial_new, period_new] = funplot_naturalPO_Hill3BP(x0, t0)
%% INPUT
% t0 : T/2
% x0 : initial
%% OUTPUT
% period_new = T/2;

%% differential correction
iteration_max = 100;
options_ODE = odeset('RelTol',1e-14, 'AbsTol',1e-14);

threshold = 1e-12;
if (abs(x0(3)) < 1e-8) && (abs(x0(6)) < 1e-8) % Lyapunov
    for s = 1 : iteration_max
        [initial_new, period_new, ~] = fun_n_differential_correction_Hill3BP_NPC_symmetry_Lyapunov(x0, t0, options_ODE);

        tspan = [0 period_new*2];
        [t_ini, x_ini] = ode113(@(t, xp) fun_ode_n_Hill3BP(t, xp), tspan, initial_new, options_ODE);
        x_ini = x_ini';

        error = abs(x_ini(:, end) - x_ini(:,1))'; % 6こ全部の誤差がthresholdを下回るように
        if error < threshold*ones(6,1)
            break;
        end

        if norm(error) > 1e+3
            disp('calculation diverged');
            return;
        end

        if s == iteration_max
            disp('do not finish');
            return;
        end
        x0 = initial_new;
        t0 = period_new; %T/2
    end
elseif (abs(x0(3)) > 1e-4) && (abs(x0(6)) < 1e-8) % halo
    for s = 1 : iteration_max
        [initial_new, period_new, ~] = fun_n_differential_correction_Hill3BP_NPC_symmetry_Halo(x0, t0, options_ODE);

        tspan = [0 period_new*2];
        [t_ini, x_ini] = ode113(@(t, xp) fun_ode_n_Hill3BP(t, xp), tspan, initial_new, options_ODE);
        x_ini = x_ini';

        error = abs(x_ini(:, end) - x_ini(:,1))'; % 6こ全部の誤差がthresholdを下回るように
        if error < threshold*ones(6,1)
            break;
        end

        if norm(error) > 1e+3
            disp('calculation diverged');
            return;
        end

        if s == iteration_max
            disp('do not finish');
            return;
        end
        x0 = initial_new;
        t0 = period_new;
    end
elseif (abs(x0(3)) < 1e-8) && (abs(x0(6)) > 1e-7) % Axial
    for s = 1 : iteration_max
        [initial_new, period_new, ~] = fun_n_differential_correction_Hill3BP_NPC_symmetry_Axial(x0, t0, options_ODE);

        tspan = [0 period_new*2];
        [t_ini, x_ini] = ode113(@(t, xp) fun_ode_n_Hill3BP(t, xp), tspan, initial_new, options_ODE);
        x_ini = x_ini';

        error = abs(x_ini(:, end) - x_ini(:,1))'; % 6こ全部の誤差がthresholdを下回るように
        if error < threshold*ones(6,1)
            break;
        end

        if norm(error) > 1e+3
            disp('calculation diverged');
            return;
        end

        if s == iteration_max
            disp('do not finish');
            return;
        end
        x0 = initial_new;
        t0 = period_new;
    end
end

%% plot
f1=figure;
axis equal
hold on
plot3(x_ini(1,:), x_ini(2,:), x_ini(3,:));
hold off;
end