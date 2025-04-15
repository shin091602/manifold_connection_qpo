function H = fun_Hamiltonian_CR3BP(xp, xt, Q, mu)
% xp  : non-dimensional position and velocity

% CR3BP
Ja = [0 1 0;
     -1 0 0;
      0 0 0];
B = [zeros(3,3); eye(3)];
R = eye(3);

r1 = sqrt((xp(1)+mu)^2+xp(2)^2+xp(3)^2);
r2 = sqrt((xp(1)-1+mu)^2+xp(2)^2+xp(3)^2);
Ux = xp(1) - (1-mu)*(xp(1)+mu)*r1^(-3) - mu*(xp(1)-1+mu)*r2^(-3);
Uy = xp(2) - (1-mu)*xp(2)*r1^(-3) - mu*xp(2)*r2^(-3);
Uz =      - (1-mu)*xp(3)*r1^(-3) - mu*xp(3)*r2^(-3);
Ur = [Ux;Uy;Uz];
f = [xp(4:6); 2*Ja*xp(4:6) + Ur];
H = (1/2)*(xp(1:6)-xt)'*Q*(xp(1:6)-xt) - (1/2)*xp(7:12)'*B*R*B'*xp(7:12) + xp(7:12)'*f; % 本当はRはR_inv、時間短くするために。
end
