function hinvcir = plot_invariant_curve(zpo,z0,p,FR)
%% plot invariatnt curve
%%% input
%zpo :periodic orbit
%z0 :initial invariant curve solution
%p :parameter

%%% output
%hinvcir :initial invariant curve plot

%% DICTIONARY OPEN
d = p("d");
N = p("N");
M = p("M");
mu = p("mu");

%% OPTIONS ODE
options = odeset('RelTol',1e-13, 'AbsTol',1e-13);

% invariant curve
X = z0(1:d*N*M);
T = z0(d*N*M+1);
rho = z0(d*N*M+2);

% fourier matrices
Fr = FR{1,2};
IFr = FR{2,2};
DFr = FR{3,2};
Q = FR{4,2};
Q = rotation_matrix(Q,rho,p);
R = IFr*Q*Fr;

% time span
inti = linspace(0,T,M+1);
inti = [inti,T];

% initilaization
YO = cell(M,N);
Yf = zeros(d*N,M);

% invariant curve calculation
for i=1:M
  % designate invariant circle
  k=d*N;
  U = X(i*k-(k-1):i*k);
  % initialize final state
  Uf = zeros(d*N,1);
  for j=1:N
    % designate initial state
    u = U(d*j-(d-1):d*j);
    % time span
    ts = [inti(i) inti(i+1)];
    % ODE
    [~,Y] = ode113(@(t,x) fun_cr3bp(t,x,mu),ts,u,options);
    % orbit state
    YO{i,j}= Y;
    % final state
    Uf(d*j-(d-1):d*j) = Y(end,1:d)';
  end
  % integrated states
  Yf(:,i) = Uf;
end

% rotation
RYf = R*Yf(:,end);

% plot invariant curve
if p("HLflg") % halo
  hinvcir = figure();
  hold on
  grid on
  axis equal
  view(10,10);
  xlabel('$x$[-]');
  ylabel('$y$[-]');
  zlabel('$z$[-]');
  % periodic orbit
  scatter3(zpo(1),zpo(2),zpo(3),400,'r.');
  % invariant circle
  for j=1:p("N")
    finvcir = scatter3(X(d*j-(d-1)),X(d*j-(d-2)),X(d*j-(d-3)),100,'m.');
    fstrobo = scatter3(Yf(d*j-(d-1),end),Yf(d*j-(d-2),end),Yf(d*j-(d-3),end),100,'b.');
    fR = scatter3(RYf(d*j-(d-1),end),RYf(d*j-(d-2),end),RYf(d*j-(d-3),end),100,'g.');
  end
  legend([finvcir,fstrobo,fR],{'$\Psi(\theta)$','$\phi_{T}(\Psi(\theta))$','$R_{-\rho}[\phi_{T}(\Psi(\theta))]$'},'Location','northoutside');
  lgd = legend;
  lgd.NumColumns = 3;
  hold off

else
  hinvcir = figure();
  hold on
  grid on
  axis equal
  xlabel('$z$[-]');
  ylabel('$\dot{z}$[-]');
  %{
  if p("congflg")
    title('GMOS, Final Iteration');
  else
    title('GMOS, First Iteration');
  end
  %}
  % periodic orbit
  scatter(zpo(3),zpo(6),400,'r.');

  % invariant circle
  for j=1:p("N")
    finvcir = scatter(X(d*j-(d-3)),X(d*j),100,'m.');
    fstrobo = scatter(Yf(d*j-(d-3),end),Yf(d*j,end),100,'b.');
    fR = scatter(RYf(d*j-(d-3),end),RYf(d*j,end),100,'g.');
  end
  legend([finvcir,fstrobo,fR],{'$\Psi(\theta)$','$\phi_{T}(\Psi(\theta))$','$R_{-\rho}[\phi_{T}(\Psi(\theta))]$'},'Location','northoutside');
  lgd = legend;
  lgd.NumColumns = 3;
  hold off
end

end