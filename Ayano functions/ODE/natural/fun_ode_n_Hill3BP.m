function  dxp = fun_ode_n_Hill3BP(t, xp)
%function [dxp, Ur, Urr] = fun_cr3bp_2n(t, xp, mu, R, Qup, Qdown)　←これはできない

% t  : non-dimensional time
% x  : non-dimensional position and velocity, x = [x y z vx vy vz]'
% xt : x_target 
% x = [r v]'
% r : non-dimensional position, r = [x y z]'
% v : non-dimensional velocity, v = [vx vy vz]'

% set 
Ja = [0 1 0;
     -1 0 0;
      0 0 0];

% xp = [r; v; pr; pv];
v = xp(4:6, :);
x = xp(1:6, :);

% the distances
r = sqrt(x(1)^2+x(2)^2+x(3)^2);

% Ur
Ur = [x(1)*(-1/r^3+3);
      -x(2)/r^3;
      -x(3)*(1/r^3+1)];

dr = v ;

dv = 2*Ja*v + Ur;

dxp = [dr;
       dv];
end