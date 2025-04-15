function Ud = Xd_finalization_qpoms(Z,Ud,p,FR)
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
mu = p("mu");

% determine finalized torus solution
% current state variables
X = Z(1:d*N*M);
T = Z(d*N*M+1);
rho = Z(d*N*M+2);

w0 = 2*pi/T;
w1 = rho/T;

% fourier matrices
Fr = FR{1,2};
IFr = FR{2,2};
DFr = FR{3,2};
Q = FR{4,2};
Q = rotation_matrix(Q,rho,p);
R = IFr*Q*Fr;

% torus function partial theta1
Ud{2,2} = DFr*X(1:d*N);

% update DUDT1 DUDT0
for i=1:N
    % location in the invariant circle
    idx = (i*d-(d-1)):i*d;
    % torus point
    U = X(idx);
    dU = fun_cr3bp([],U,mu);
    Ud{1,2}(idx) = (1/w0).*(dU-w1*Ud{2,2}(idx));
end

% update previous solutions
Ud{3,2} = Z;

end