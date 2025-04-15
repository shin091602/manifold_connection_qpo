function [A_11, V_n, D_n, eig_n, real_D_n, imag_D_n] = fun_n_A_Hill3BP(r0)
% 平衡点の座標を与えると、状態行列Aと固有値、固有ベクトルを返してくれる関数
% hill3BP
% Quadratic optimal control. パラメーターはQ(R=I).

% hill3BP
Ja = [0 1 0;
     -1 0 0;
      0 0 0];
r = sqrt(r0(1)^2+r0(2)^2+r0(3)^2);
   
%% calculate A
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

A_11 = [zeros(3) eye(3);
        Urr 2*Ja];

[V_n,D_n] = eig(A_11); % V=固有ベクトル, D=固有値

eig_n = zeros(1, 6);
for i = 1 : 6
    eig_n(i) = D_n(i,i);
end

real_D_n = zeros(1,6);
imag_D_n = zeros(1,6);

for i = 1:6
    real_D_n(i) = real(D_n(i,i));
    imag_D_n(i) = imag(D_n(i,i));
end

