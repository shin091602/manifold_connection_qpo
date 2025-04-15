function [A, V_n, D_n, eig_n, real_D_n, imag_D_n] = fun_n_A_CR3BP(r0, mu)
% 平衡点の座標とmuを与えると、状態行列Aと固有値、固有ベクトルを返してくれる関数
% CR3BP

Ja = [0 1 0;
     -1 0 0;
      0 0 0];
   
%% calculate A
r1 = sqrt((r0(1)+mu)^2+r0(2)^2+r0(3)^2);
r2 = sqrt((r0(1)-1+mu)^2+r0(2)^2+r0(3)^2);

Ux = r0(1) - (1-mu)*(r0(1)+mu)*r1^(-3) - mu*(r0(1)-1+mu)*r2^(-3);
Uy = r0(2) - (1-mu)*r0(2)*r1^(-3) - mu*r0(2)*r2^(-3);
Uz =       - (1-mu)*r0(3)*r1^(-3) - mu*r0(3)*r2^(-3);

Ur = [Ux;Uy;Uz];

Uxx = 1-(1-mu)*(r1^(-3)-3*(r0(1)+mu)^2*r1^(-5)) - mu*(r2^(-3)-3*(r0(1)-1+mu)^2*r2^(-5));
Uxy = 3*r0(2)*((1-mu)*(r0(1)+mu)*r1^(-5) + mu*(r0(1)-1+mu)*r2^(-5));
Uxz = 3*r0(3)*((1-mu)*(r0(1)+mu)*r1^(-5) + mu*(r0(1)-1+mu)*r2^(-5));
Uyy = 1 - (1-mu)*(r1^(-3)-3*r0(2)^2*r1^(-5)) - mu*(r2^(-3) - 3*r0(2)^2*r2^(-5));
Uyz = 3*r0(2)*r0(3)*((1-mu)*r1^(-5) + mu*r2^(-5));
Uzz = -(1-mu)*(r1^(-3) - 3*r0(3)^2*r1^(-5)) - mu*(r2^(-3) - 3*r0(3)^2*r2^(-5));

Urr = [Uxx Uxy Uxz;
       Uxy Uyy Uyz;
       Uxz Uyz Uzz];

A = [zeros(3) eye(3);
        Urr 2*Ja];

[V_n,D_n] = eig(A); % V=固有ベクトル, D=固有値

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

