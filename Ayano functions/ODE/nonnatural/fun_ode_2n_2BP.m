function  dxp = fun_ode_2n_2BP(t, xp, Q, xt, mu)
% t  : non-dimensional time
% x  : non-dimensional position and velocity, x = [x y z vx vy vz]'
% xt : x_target 
% x = [r v]'
% r : non-dimensional position, r = [x y z]'
% v : non-dimensional velocity, v = [vx vy vz]'
% p : p = [pr pv]'
% pr = [px py pz]'
% pv = [pxdot pydot pzdot]'

% set
Ja = zeros(3,3);
R = eye(3);
Qup = Q(1:3, :);
Qdown = Q(4:6, :);

% xp = [r; v; pr; pv];
v = xp(4:6, :);
pr = xp(7:9, :);
pv = xp(10:12, :);
x = xp(1:6, :);
r0 = xp(1:3, :);

% the distances
r = sqrt(x(1)^2+x(2)^2+x(3)^2);

% Ur ok
Ur = -mu/r^3*[r0(1); r0(2); r0(3)];

% Urr ok
rTr = [r0(1)^2 r0(1)*r0(2) r0(1)*r0(3);
       r0(1)*r0(2) r0(2)^2 r0(2)*r0(3);
       r0(1)*r0(3) r0(2)*r0(3) r0(3)^2];

Urr = mu/r^5*(3*rTr-r^2*eye(3));

R_inv = inv(R);

dr = v ;

dv = 2*Ja*v + Ur - R_inv*pv;

dpr = -Qup*(x-xt) - Urr*pv;

dpv = -Qdown*(x-xt) - pr + 2*Ja*pv;

dxp = [dr;
       dv;
       dpr;
       dpv];
end