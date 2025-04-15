function  dxp = fun_ode_2n_CR3BP(t, xp, mu, Q, xt)
%function [dxp, Ur, Urr] = fun_cr3bp_2n(t, xp, mu, R, Qup, Qdown)　←これはできない

% t  : non-dimensional time
% x  : non-dimensional position and velocity, x = [x y z vx vy vz]'
% xt : x_target 
% x = [r v]'
% r : non-dimensional position, r = [x y z]'
% v : non-dimensional velocity, v = [vx vy vz]'
% p : p = [pr pv]'
% pr = [p1 p2 p3]'
% pv = [p1dot p2dot p3dot]'
% mu : mass ratio of the primaries

% set
Ja = [0 1 0;
     -1 0 0;
      0 0 0];
R = eye(3);
Qup = Q(1:3, :);
Qdown = Q(4:6, :);

% xp = [r; v; pr; pv];
r = xp(1:3, :);
v = xp(4:6, :);
pr = xp(7:9, :);
pv = xp(10:12, :);
x = xp(1:6, :);

%% calculate A
r1 = sqrt((r(1)+mu)^2+r(2)^2+r(3)^2);
r2 = sqrt((r(1)-1+mu)^2+r(2)^2+r(3)^2);

Ux = r(1) - (1-mu)*(r(1)+mu)*r1^(-3) - mu*(r(1)-1+mu)*r2^(-3);
Uy = r(2) - (1-mu)*r(2)*r1^(-3) - mu*r(2)*r2^(-3);
Uz =      - (1-mu)*r(3)*r1^(-3) - mu*r(3)*r2^(-3);

Ur = [Ux;Uy;Uz];

Uxx = 1-(1-mu)*(r1^(-3)-3*(r(1)+mu)^2*r1^(-5)) - mu*(r2^(-3)-3*(r(1)-1+mu)^2*r2^(-5));
Uxy = 3*r(2)*((1-mu)*(r(1)+mu)*r1^(-5) + mu*(r(1)-1+mu)*r2^(-5));
Uxz = 3*r(3)*((1-mu)*(r(1)+mu)*r1^(-5) + mu*(r(1)-1+mu)*r2^(-5));
Uyy = 1 - (1-mu)*(r1^(-3)-3*r(2)^2*r1^(-5)) - mu*(r2^(-3) - 3*r(2)^2*r2^(-5));
Uyz = 3*r(2)*r(3)*((1-mu)*r1^(-5) + mu*r2^(-5));
Uzz = -(1-mu)*(r1^(-3) - 3*r(3)^2*r1^(-5)) - mu*(r2^(-3) - 3*r(3)^2*r2^(-5));

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