function [F_poms,Xd] = F_poms(z_mul,Xd,p)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Error vector for multiple shooting constraints and phase conditions
%% By:Soi Yamaguchi
%%% input
% z_mul(1:d*M) :X of patch points
% z_mul(d*M+1) :the period of the periodic orbit
% Xd :data package :x0,f0
% p:paramter dictionary

%%% output
% F_poms:Error vector
% Xd: updated data package:x0,f0,xT,phiT,f
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% DICTIONARY
d = p("d");
M = p("M");
mu = p("mu");
Cref = p("C");

%% OPTIONS ODE
options = odeset('RelTol',1e-13, 'AbsTol',1e-13);

% period of periodic orbit
T = z_mul(d*M+1);
% tspan
tspan = linspace(0, T-1e-5, M+1);
% initialization
xT = zeros(d,M);
phiT = zeros(d,d,M);
f = zeros(d,M);

% patch points for multiple shooting
for i=1:M
  if i==M
    tspan_pt = [tspan(i) T];
  else
    tspan_pt = [tspan(i) tspan(i+1)];
  end
  X0 = [z_mul(d*i-(d-1):d*i); reshape(eye(d), d^2, 1)];
  [~,Y] = ode113(@(t,x) fun_stm_cr3bp(t,x,mu),tspan_pt,X0,options);
  % final state
  xT(:,i) = Y(end,1:d)';
  % final STM
  phiT(:,:,i) = reshape(Y(end,d+1:d+d^2),d,d);
  % vector field at final state
  f(:,i) = fun_cr3bp([],xT(:,i),mu);
end

% add to data package
Xd{3,1} = "xT";
Xd{3,2} = xT;
Xd{4,1} = "phiT";
Xd{4,2} = phiT;
Xd{5,1} = "f";
Xd{5,2} = f;

% error vector
F = zeros(d*M,1);
for i=2:M
  j=i-1;
  F(d*j-(d-1):d*j) = xT(:,j)-z_mul(d*i-(d-1):d*i);
end
F((M-1)*d+1:M*d) = xT(:,end)-z_mul(1:d);

% parametrization constraints
Ccon_po = Jacobi_const(z_mul(1:d),mu)-Cref;

% phase condition constraints
x0 = Xd{1,2};
f0 = Xd{2,2};
phs = dot(z_mul(1:d)-x0,f0);

% F_poms
F_poms = [F;phs;Ccon_po];

end