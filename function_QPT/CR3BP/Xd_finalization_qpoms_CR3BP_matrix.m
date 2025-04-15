function Ud = Xd_finalization_qpoms_CR3BP_matrix(Z,Ud,p)
%%% input
% Zd :torus solution
% Ud: finalized data package
% p:variables in dictionary

%%% output
% Xd: updated data package:xT,phiT,f,x0,f0

%% DICTIONARY OPEN
d = p("d");
N = p("N");
M = p("M");

% determine finalized torus solution
% current state variables
X = Z(1:d*N*M);
T = Z(d*N*M+1);
rho = Z(d*N*M+2);

w0 = 2*pi/T;
w1 = rho/T;

% fourier matrices
[~, DFR_theta, ~] = fun_Fourier_matrix(rho,p);

% torus function partial theta1
Ud{2,2} = DFR_theta*X(1:d*N);

% update DUDT1 DUDT0
for i=1:N
    % location in the invariant circle
    idx = (i*d-(d-1)):i*d;
    % torus point
    U = X(idx);
    dU = fun_ode_n_CR3BP([],U,p("mu"));
    Ud{1,2}(idx) = (1/w0).*(dU-w1*Ud{2,2}(idx));
end

% update previous solutions
Ud{3,2} = Z;

end