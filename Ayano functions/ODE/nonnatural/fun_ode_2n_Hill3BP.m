function  dxp = fun_ode_2n_Hill3BP(t, xp, Q, xt)

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
Ja = [0 1 0;
     -1 0 0;
      0 0 0];
R = eye(3);
Qup = Q(1:3, :);
Qdown = Q(4:6, :);

% xp = [r; v; pr; pv];
v = xp(4:6, :);
pr = xp(7:9, :);
pv = xp(10:12, :);
x = xp(1:6, :);

% the distances
r = sqrt(x(1)^2+x(2)^2+x(3)^2);

% Ur
Ur = [x(1)*(-1/r^3+3);
      -x(2)/r^3;
      -x(3)*(1/r^3+1)];

% Urr
Uxx = 3 - 1/r^3 + 3*x(1)^2/r^5;
Uxy = 3*x(1)*x(2)/r^5;
Uxz = 3*x(1)*x(3)/r^5;
Uyy = -1/r^3 + 3*x(2)^2/r^5;
Uyz = 3*x(2)*x(3)/r^5;
Uzz = -1 -1/r^3 + 3*x(3)^2/r^5;

Urr = [Uxx Uxy Uxz;
       Uxy Uyy Uyz; 
       Uxz Uyz Uzz];

dr = v ;

dv = 2*Ja*v + Ur - R\pv;

dpr = -Qup*(x-xt) - Urr*pv;

dpv = -Qdown*(x-xt) - pr + 2*Ja*pv;

dxp = [dr;
       dv;
       dpr;
       dpv];
end