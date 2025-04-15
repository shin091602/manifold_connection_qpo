function [A, V_2n, D_2n, eig_2n, real_D_2n, imag_D_2n, Ur, Urr, A_11] = fun_2n_A_2BP(Q, r0, pv0, mu)
% 平衡点の座標とQを与えると、状態行列Aと固有値、固有ベクトルを返してくれる関数
% 2BP
% Quadratic optimal control. パラメーターはQ(R=I).

B = [zeros(3); eye(3)];
R = eye(3)*10^(0);
R_inv = inv(R);
Q_1313 = Q(1:3, 1:3);
Q_1346 = Q(1:3, 4:6);
Qup = Q(1:3, 1:6);
Qdown = Q(4:6, 1:6);

% 2BP
Ja = [0 0 0;
      0 0 0;
      0 0 0];
   
%% calculate A
% r
r = sqrt(r0(1)^2+r0(2)^2+r0(3)^2);

% Ur ok
Ur = -mu/r^3*[r0(1); r0(2); r0(3)];

% Urr ok
rTr = [r0(1)^2 r0(1)*r0(2) r0(1)*r0(3);
       r0(1)*r0(2) r0(2)^2 r0(2)*r0(3);
       r0(1)*r0(3) r0(2)*r0(3) r0(3)^2];

Urr = mu/r^5*(3*rTr-r^2*eye(3));

% Urrx ok
Uxxx = 3*mu*r0(1)/r^7*(-2*r0(1)^2+3*r0(2)^2+3*r0(3)^2); 
Uxyx = 3*mu*r0(2)/r^7*(-4*r0(1)^2+r0(2)^2+r0(3)^2); 
Uxzx = 3*mu*r0(3)/r^7*(-4*r0(1)^2+r0(2)^2+r0(3)^2);
Uyyx = 3*mu*r0(1)/r^7*(r0(1)^2-4*r0(2)^2+r0(3)^2);
Uyzx = -15*mu*r0(1)*r0(2)*r0(3)/r^7;
Uzzx = 3*mu*r0(1)/r^7*(r0(1)^2+r0(2)^2-4*r0(3)^2);

Urrx=[Uxxx Uxyx Uxzx;
      Uxyx Uyyx Uyzx;
      Uxzx Uyzx Uzzx];

% Urry ok
Uxxy = 3*mu*r0(2)/r^7*(-4*r0(1)^2+r0(2)^2+r0(3)^2); 
Uxyy = 3*mu*r0(1)/r^7*(r0(1)^2-4*r0(2)^2+r0(3)^2); 
Uxzy = -15*mu*r0(1)*r0(2)*r0(3)/r^7;
Uyyy = 3*mu*r0(2)/r^7*(3*r0(1)^2-2*r0(2)^2+3*r0(3)^2);
Uyzy = 3*mu*r0(3)/r^7*(r0(1)^2-4*r0(2)^2+r0(3)^2);
Uzzy = 3*mu*r0(2)/r^7*(r0(1)^2+r0(2)^2-4*r0(3)^2);

Urry=[Uxxy Uxyy Uxzy;
      Uxyy Uyyy Uyzy;
      Uxzy Uyzy Uzzy];

% Urrz ok
Uxxz = 3*mu*r0(3)/r^7*(-4*r0(1)^2+r0(2)^2+r0(3)^2);
Uxyz = -15*mu*r0(1)*r0(2)*r0(3)/r^7; 
Uxzz = 3*mu*r0(1)/r^7*(r0(1)^2+r0(2)^2-4*r0(3)^2);
Uyyz = 3*mu*r0(3)/r^7*(r0(1)^2-4*r0(2)^2+r0(3)^2);
Uyzz = 3*mu*r0(2)/r^7*(r0(1)^2+r0(2)^2-4*r0(3)^2);
Uzzz = 3*mu*r0(3)/r^7*(3*r0(1)^2+3*r0(2)^2-2*r0(3)^2);

Urrz=[Uxxz Uxyz Uxzz;
      Uxyz Uyyz Uyzz;
      Uxzz Uyzz Uzzz];

Urrrpv0=[Urrx*pv0 Urry*pv0 Urrz*pv0];

A_11 = [zeros(3) eye(3);
        Urr 2*Ja];
A_12 = [zeros(3) zeros(3);
        zeros(3) -R_inv];
A_21 = [-Q_1313-Urrrpv0 -Q_1346;
        -Qdown];
A_22 = [zeros(3) -Urr;
        -eye(3) 2*Ja];

% A_21 = [-Q_1313 -Q_1346;
%         -Qdown];

A = [A_11 A_12;
     A_21 A_22];

[V_2n,D_2n] = eig(vpa(A)); % V=固有ベクトル, D=固有値

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

% [V_2n,D_2n] = eig(vpa(A)); % V=固有ベクトル, D=固有値 % vpa : 可変精度の演算を使用して 5 次の魔方陣行列の数値としての固有値を計算します。



