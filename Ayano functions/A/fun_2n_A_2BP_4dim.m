function [A, V_2n, D_2n, eig_2n, real_D_2n, imag_D_2n, Ur, Urr] = fun_2n_A_2BP_4dim(Q, r0, pv0, mu)
% 平衡点の座標とQを与えると、状態行列Aと固有値、固有ベクトルを返してくれる関数
% 2BP
% Quadratic optimal control. パラメーターはQ(R=I).

% set
R_inv = 1;
Ja = 0;
   
%% calculate A
% r
r = sqrt(r0(1)^2+r0(2)^2+r0(3)^2);

% Ur ok
Ur = -mu/r^2;

% Urr ok
Urr = 2*mu/r^3;

% Urrr ok
Urrr=-6*mu/r^4;

A_11 = [0 1;
        Urr 2*Ja];
A_12 = [0 0;
        0 -R_inv];
A_21 = [-Q(1,1)-Urrr*pv0 -Q(1,2);
        -Q(2,1) -Q(2,2)];
A_22 = [0 -Urr;
        -1 2*Ja];

A = [A_11 A_12;
     A_21 A_22];

[V_2n,D_2n] = eig(vpa(A)); % V=固有ベクトル, D=固有値

eig_2n = zeros(1, 4);
for i = 1 : 4
    eig_2n(i) = D_2n(i,i);
end

real_D_2n = zeros(1,4);
imag_D_2n = zeros(1,4);

for i = 1:4
    real_D_2n(i) = real(D_2n(i,i));
    imag_D_2n(i) = imag(D_2n(i,i));
end

% [V_2n,D_2n] = eig(vpa(A)); % V=固有ベクトル, D=固有値 % vpa : 可変精度の演算を使用して 5 次の魔方陣行列の数値としての固有値を計算します。



