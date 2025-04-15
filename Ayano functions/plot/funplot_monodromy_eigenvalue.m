function [fig_monodromy_eig, fig_monodromy_eig_zoom] = funplot_monodromy_eigenvalue(eig, vec)

%% plot modromy's eigenvalue
fig_monodromy_eig = figure;
hold on;
axis equal;

theta = linspace(0, 2*pi, 1001)';
x_circle = cos(theta);
y_circle = sin(theta);
plot(x_circle,y_circle,'k');
% plot orbit
for i=1:length(eig)
    if (abs(vec(3, i)) < 1e-8) && (abs(vec(6, i)) < 1e-8) % inplane
        plot(real(eig(i)), imag(eig(i)), '*', 'color', 'b', 'MarkerSize',10);
    elseif (abs(vec(1, i)) < 1e-8) && (abs(vec(2, i)) < 1e-8) && (abs(vec(4, i)) < 1e-8) && (abs(vec(5, i)) < 1e-8) % out-plane
        plot(real(eig(i)), imag(eig(i)), '*', 'color', 'r', 'MarkerSize',10);
    else % in-plane + out-plane
        plot(real(eig(i)), imag(eig(i)), '*', 'color', 'g', 'MarkerSize',10);
    end
end
xlabel('Re[$\lambda_i$]');
ylabel('Im[$\lambda_i$]');
grid on;
hold off;

%% plot modromy's eigenvalue (zoom)
fig_monodromy_eig_zoom = figure;
hold on;
axis equal;

plot(x_circle,y_circle,'k');
% plot orbit
for i=1:length(eig)
    if (abs(vec(3, i)) < 1e-8) && (abs(vec(6, i)) < 1e-8) % inplane
        plot(real(eig(i)), imag(eig(i)), '*', 'color', 'b', 'MarkerSize',10);
    elseif (abs(vec(1, i)) < 1e-8) && (abs(vec(2, i)) < 1e-8) && (abs(vec(4, i)) < 1e-8) && (abs(vec(5, i)) < 1e-8) % out-plane
        plot(real(eig(i)), imag(eig(i)), '*', 'color', 'r', 'MarkerSize',10);
    else % in-plane + out-plane
        plot(real(eig(i)), imag(eig(i)), '*', 'color', 'g', 'MarkerSize',10);
    end
end
xlabel('Re[$\lambda_i$]');
ylabel('Im[$\lambda_i$]');
xlim([-1.5 1.5]);
ylim([-1.5 1.5]);
grid on;
hold off;
end