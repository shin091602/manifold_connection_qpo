function Q = rotation_matrix(Q,rho,p)
% fourier transformation
%%% input
%rho :rotation number
%d :dimension
%N :the number of points in invariant circle

%%% output
%Q :rotation matrix

%% DICTIONARY OPEN
d = p("d");
N = p("N")+1;

% order -remove K=0
if mod(N,2)~=0
    K = -(N-1)/2:1:(N-1)/2;
    for i=1:length(K)
        if K(i)==0
            J = i;
        end
    end
    K = [K(J:end),K(1:J-1)]';
else
    K=-N/2:1:N/2;
    for i=1:length(K)
        if K(i)==0
            J = i;
        end
    end
    K = [K(J:end),K(1:J-1)]';
end

% initialize
k = 1;
for i=2:2:(N-1)
    %column index
    k = k+1;
    %indexing
    row1 = d*i-(d-1):d*i;
    row2 = d*i+1:d*(i+1);
    %rotating fourier coefficients
    Q(row1,row1) = cos(K(k)*rho).*eye(d);
    Q(row1,row2) = -sin(K(k)*rho).*eye(d);
    Q(row2,row1) = sin(K(k)*rho).*eye(d);
    Q(row2,row2) = cos(K(k)*rho).*eye(d);
end
% Add first block
Q(1:d, 1:d) = eye(d);


end