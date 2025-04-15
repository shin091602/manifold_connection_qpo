function [z_mul,Xd] = po_mulshoot_initialization(z0,p)
%%% multiple shooting patch points at the periodic orbit
%%% input
% z0(1:d) :X of patch points
% z0(d+1) :the period of the periodic orbit
% p:variables in dictionary

%%% output
% z_mul(1:d*M):X of patch points
% Xd:data package:fT,x0,f0

%% OPTIONS ODE
options = odeset('RelTol',1e-13, 'AbsTol',1e-13);

% dictionary 
d = p("d");
M = p("M");
mu = p("mu");

% patch times
pt = linspace(0,z0(d+1),M+1);
% solve ODE
[~,x] = ode113(@(t,x) fun_cr3bp(t,x,mu),pt,z0(1:d),options);

% arrange patch points' data
z_mul = zeros(d*M+1,1);
for i=1:length(pt)-1
z_mul(d*i-5:d*i) = x(i,:);
end
%orbital period
z_mul(d*M+1) = z0(d+1);

% data package
data_size = 5;
Xd = cell(data_size,2);%[name,values]
Xd{1,1} = "x0";
Xd{1,2} = z_mul(1:d);
Xd{2,1} = "f0";
Xd{2,2} = fun_cr3bp([],z_mul(1:d),mu);

end