function [hPO,X] = plot_periodic_orbit(z,p,L2)
%%% plot figure of refined periodic %%% input
%z:orbita data of periodic orbit
%p:variable dictionary

%%% output
%hPO :figure of two orbits(refined PO--blue, initial guessed PO--black)
%X :total orbit data

%options ODE
options_ODE = odeset('RelTol',1e-13, 'AbsTol',1e-13);

%figure of periodic orbit
hPO = figure();
hold on
grid on
box on
axis equal
xlabel('$x$[-]');
ylabel('$y$[-]');
zlabel('$z$[-]');
view([-22 21]);
%view([0 90]);
%libration point
plot3(L2(1),L2(2),L2(3),'*','MarkerFaceColor','black','MarkerEdgeColor','black','MarkerSize',8);

% tspan
% tspan = linspace(0, z(end)-1e-5, p("M")+1);
tspan = linspace(0, z(end), p("M")+1);
for i=1:p("M")
  if i==p("M")
    tspan_pt = [tspan(i) z(end)];
  else
    tspan_pt = [tspan(i) tspan(i+1)];
  end
  [~,X] = ode113(@(t,x) fun_cr3bp(t,x,p("mu")),tspan_pt,z(p("d")*i-5:p("d")*i)',options_ODE);
  %patch points plot
  if i~=p("M")
  fpatch = scatter3(z(p("d")*i-5),z(p("d")*i-4),z(p("d")*i-3),50,'bo','filled');
  else
    fpatch = scatter3(z(p("d")*i-5),z(p("d")*i-4),z(p("d")*i-3),50,'ro','filled');
  end

  %orbit plot
  fperi = plot3(X(:,1),X(:,2),X(:,3),'k','LineWidth',1.5);
  hold on
end
lgd = legend([fperi,fpatch],{'Periodic orbit','Patch points'},'Location','northeast');
lgd.NumColumns=1;
hold off

end