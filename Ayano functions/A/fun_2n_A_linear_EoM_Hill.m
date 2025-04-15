function [A, V_2n, D_2n, eig_2n, real_D_2n, imag_D_2n] = fun_2n_A_linear_EoM_Hill(Q, r0)
% 平衡点の座標とQを与えると、状態行列Aと固有値、固有ベクトルを返してくれる関数
% hill3BP
% Quadratic optimal control. パラメーターはQ(R=I).
%% Urrrpvの項がないとこうなる．

B = [zeros(3); eye(3)];
R = eye(3)*10^(0);
R_inv = inv(R);
Q_1313 = Q(1:3, 1:3);
Q_1346 = Q(1:3, 4:6);
Qup = Q(1:3, 1:6);
Qdown = Q(4:6, 1:6);

% hill3BP
Ja = [0 1 0;
     -1 0 0;
      0 0 0];
   
%% calculate A
r = sqrt(r0(1)^2+r0(2)^2++r0(3)^2);

% Ur
Ur = [r0(1)*(-1/r^3+3);
      -r0(2)/r^3;
      -r0(3)*(1/r^3+1)];

% Urr
Uxx = 3 - 1/r^3 + 3*r0(1)^2/r^5;
Uxy = 3*r0(1)*r0(2)/r^5;
Uxz = 3*r0(1)*r0(3)/r^5;
Uyy = -1/r^3 + 3*r0(2)^2/r^5;
Uyz = 3*r0(2)*r0(3)/r^5;
Uzz = -1 -1/r^3 + 3*r0(3)^2/r^5;

Urr = [Uxx Uxy Uxz;
       Uxy Uyy Uyz; 
       Uxz Uyz Uzz];
% % Urrx
% Uxxx = 3*r0(1)/r^5*(3-5*r0(1)^2/r^2);
% Uxyx = 3*r0(2)/r^5*(1-5*r0(1)^2/r^2);
% Uxzx = 3*r0(3)/r^5*(1-5*r0(1)^2/r^2);
% Uyyx = 3*r0(1)/r^5*(1-5*r0(2)^2/r^2);
% Uyzx = -15*r0(1)*r0(2)*r0(3)/r^7;
% Uzzx = 3*r0(1)/r^5*(1-5*r0(3)^2/r^2);
% 
% Urrx=[Uxxx Uxyx Uxzx;
%       Uxyx Uyyx Uyzx;
%       Uxzx Uyzx Uzzx];
% 
% % Urry
% Uxxy = 3*r0(2)/r^5*(1-5*r0(1)^2/r^2);
% Uxyy = 3*r0(1)/r^5*(1-5*r0(2)^2/r^2);
% Uxzy = -15*r0(1)*r0(2)*r0(3)/r^7;
% Uyyy = 3*r0(2)/r^5*(3-5*r0(2)^2/r^2);
% Uyzy = 3*r0(3)/r^5*(1-5*r0(2)^2/r^2);
% Uzzy = 3*r0(2)/r^5*(1-5*r0(3)^2/r^2);
% 
% Urry=[Uxxy Uxyy Uxzy;
%       Uxyy Uyyy Uyzy;
%       Uxzy Uyzy Uzzy];
% 
% % Urrz
% Uxxz = 3*r0(3)/r^5*(1-5*r0(1)^2/r^2);
% Uxyz = -15*r0(1)*r0(2)*r0(3)/r^7;
% Uxzz = 3*r0(1)/r^5*(1-5*r0(3)^2/r^2);
% Uyyz = 3*r0(3)/r^5*(1-5*r0(2)^2/r^2);
% Uyzz = 3*r0(2)/r^5*(1-5*r0(3)^2/r^2);
% Uzzz = 3*r0(3)/r^5*(3-5*r0(3)^2/r^2);
% 
% Urrz=[Uxxz Uxyz Uxzz;
%       Uxyz Uyyz Uyzz;
%       Uxzz Uyzz Uzzz];

% Urrrpv0=[Urrx*pv0 Urry*pv0 Urrz*pv0];

A_11 = [zeros(3) eye(3);
        Urr 2*Ja];
A_12 = [zeros(3) zeros(3);
        zeros(3) -R_inv];
A_21 = [-Q_1313 -Q_1346;
        -Qdown];
A_22 = [zeros(3) -Urr;
        -eye(3) 2*Ja];

% A_21 = [-Q_1313 -Q_1346;
%         -Qdown];

A = [A_11 A_12;
     A_21 A_22];

[V_2n,D_2n] = eig(A); % V=固有ベクトル, D=固有値
% [V_n, D_n] = eig(A_11);

eig_2n = zeros(1, 12);
for i = 1 : 12
    eig_2n(i) = D_2n(i,i);
end

real_D_2n = zeros(1,12);
imag_D_2n = zeros(1,12);

for i = 1:12
    real_D_2n(i) = real(D_2n(i,i));
    imag_D_2n(i) = imag(D_2n(i,i));
end

