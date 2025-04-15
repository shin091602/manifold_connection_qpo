function [R,IR,DR] = fourier_matrix(p)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% fourier transformation
%% By:Soi Yamaguchi
%%% input
%d :dimension
%N :the number of points in invariant circle

%%% output
%R :fourier trans
%IR :inverse fourier trans
%DR :derivative matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% DICTIONARY OPEN
d = p("d");
N = p("N");

% initialize % matrix_size=90*90
D = zeros(d*N,d*N);R = zeros(d*N,d*N);
IR = zeros(d*N,d*N);DR = zeros(d*N,d*N);

% order -remove K=0
K = -(N-1)/2:1:(N-1)/2; % -7:1:7(15点)
for i=1:length(K)
    if K(i)==0
        J = i;
    end
end
K = [K(J:end),K(1:J-1)]'; % K=[0 1 2 3 4 5 6 7 -7 -6 -5 -4 -3 -2 -1]なぜ位置地並び替えたかは不明

for i=1:N
    %k index
    k =1;
    %row index
    rdx = d*(i-1)+1:d*i;
    for j=1:N
        %column index
        cdx = d*(j-1)+1:d*j;
        D(rdx,cdx) = eye(d).*exp(2*pi*K(i)*K(j)*(1i)/N);
        %inverse fourier matrix
        if mod(j,2) ==0
            k = k+1;
            IR(rdx,cdx) = cos(2*pi*(i-1)*K(k)/N)*eye(d);
            DR(rdx,cdx) = -K(k)*sin(2*pi*(i-1)*K(k)/N)*eye(d);
        else
            IR(rdx,cdx) = sin(2*pi*(i-1)*K(k)/N)*eye(d);
            DR(rdx,cdx) = K(k)*cos(2*pi*(i-1)*K(k)/N)*eye(d);
        end
    end
    IR(rdx,1:d) = eye(d);
end

% complex to real
k=0;
for i=2:2:N-1
    %column index
    k = k+1;
    %indexing
    r1 = (d*(i-1)+1):d*i;
    r2 = (d*i+1):d*(i+1);
    c1 = (d*k+1):d*(k+1);
    c2 = (d*(N-k)+1):d*(N-k+1);
    %define R
    R(r1,c1) = eye(d)/N;
    R(r1,c2) = eye(d)/N;
    R(r2,c1) = -1i.*eye(d)/N;
    R(r2,c2) = 1i.*eye(d)/N;
end
R(1:d,1:d) = eye(d)/N;
R = R*D;
if abs(max(imag(R)))>1e-16
    disp("R is not exactly real")
end
R = real(R);

% derivative matrix
DR = DR*R;
if abs(max(imag(R)))>1e-16
    disp("R is not exactly real")
end
DR = real(DR);

end