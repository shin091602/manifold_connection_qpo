function  dxp = fun_ode_4dim_2BP(t, xp, Q, xt, mu)
%　2BPの4次元のとき

% t  : non-dimensional time
% x  : non-dimensional position and velocity, x = [x xdot]'
% xt : x_target 
% r : non-dimensional position, r = [x]'
% v : non-dimensional velocity, v = [xdot]'
% p : p = [pr pv]'
% pr = [px]'
% pv = [pxdot]'

% set
Ja = [0 1 0;
     -1 0 0;
      0 0 0];
R = eye(3);
Qup = Q(1:3, :);
Qdown = Q(4:6, :);

% xp = [x; xdot; px; pxdot];
v = xp(2, :);
pr = xp(3, :);
pv = xp(4, :);
x = xp(1:2, :);
r0 = xp(1, :);

% the distances
r = xp(1);

% Ur ok
Ur = -mu/xp(1)^2;

% Urr ok
Urr = 2*mu/xp(1)^3;

R_inv = inv(R);

dr = v;

dv = 2*Ja*v + Ur - R_inv*pv;

dpr = -Qup*(x-xt) - Urr*pv;

dpv = -Qdown*(x-xt) - pr + 2*Ja*pv;

dxp = [dr;
       dv;
       dpr;
       dpv];
end