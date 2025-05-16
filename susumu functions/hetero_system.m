function F = hetero_system(th, interpLyap1, interpLyap2, vU1, vS2, p, tspan_u, tspan_s, opt1, opt2)
% th(1): θ₁ (L1 軌道上の角度)
% th(2): θ₂ (L2 軌道上の角度)

% --- L1 の不安定多様体 ---
X1     = interpLyap1(th(1));
x1     = X1(1:6);
phi1   = reshape(X1(7:end), 6, 6);
U      = phi1 * vU1;
U = U/norm(U);
XU     = x1' + p('pert') * U;
[~,Y1,~,YE1,~] = ode113(@(t,x) fun_cr3bp(t,x,p('mu')), tspan_u, XU, opt1);
if ~isempty(YE1)
    P1 = YE1(end, :);
else
    P1 = Y1(end, :);
end

% --- L2 の安定多様体（逆方向伝播）---
X2     = interpLyap2(th(2));
x2     = X2(1:6);
phi2   = reshape(X2(7:end), 6, 6);
S      = phi2 * vS2;
S = S/norm(S);
XS     = x2' + p('pert') * S;
[~,Y2,~, YE2,~] = ode113(@(t,x) fun_cr3bp(t,x,p('mu')), tspan_s, XS, opt2);
if ~isempty(YE2)
    P2 = YE2(end, :);
else
    P2 = Y2(end, :);
end

% ２つの断面交点差 (2×1 ベクトル)
F = (P1([2, 4]) - P2([2, 4]))';
end