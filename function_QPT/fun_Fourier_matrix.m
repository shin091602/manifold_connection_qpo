function [FR, DFR_theta, DFR_rho] = fun_Fourier_matrix(rho, p)

%% DICTIONARY OPEN
d = p("d");
N = p("N");

%% make DFT matrix
ck = zeros(N,N);
for k = 1:N % k:周波数, row
    K = -(N-1)/2-1+k;
    for j = 1:N % j=j+t/T, j=0~N-1, column
        ck(k,j) = exp(-2*pi*(1i)*K*(j-1)/N); %eq(3.38)
    end
end
DFT = kron(ck, eye(d)); %離散フーリエ変換

%% make IDFT matrix
Ick_Rho = zeros(N,N);
for j = 1:N % row
    for k = 1:N % column
        K = -(N-1)/2-1+k;
        Ick_Rho(j,k) = 1/N*exp(2*pi*(1i)*K*(j-1)/N)*exp(-1i*K*rho);% 逆離散フーリエ変換
    end
end
IDFT_rho = kron(Ick_Rho, eye(d)); %離散フーリエ変換

Ick = zeros(N,N);
for j = 1:N % row
    for k = 1:N % column
        K = -(N-1)/2-1+k;
        Ick(j,k) = 1/N*exp(2*pi*(1i)*K*(j-1)/N);% 逆離散フーリエ変換
    end
end
check_DFT = Ick*ck;
for i = 1:N
    for j = 1:N
        if abs(imag(check_DFT(i,j))) < 1e-13
            check_DFT(i,j) = real(check_DFT(i,j));
        end
        if abs(check_DFT(i,j)) < 1e-13
            check_DFT(i,j) = 0;
        end
    end
end

%% make FR
FR = IDFT_rho*DFT;
for i = 1:d*N
    for j = 1:d*N
        if abs(imag(FR(i,j))) < 1e-13
            FR(i,j) = real(FR(i,j));
        end
        if abs(FR(i,j)) < 1e-13
            FR(i,j) = 0;
        end
    end
end

%% make FR_check
IDFT = kron(Ick, eye(d)); %離散フーリエ変換
FR_check = IDFT*DFT;
for i = 1:d*N
    for j = 1:d*N
        if abs(imag(FR_check(i,j))) < 1e-13
            FR_check(i,j) = real(FR_check(i,j));
        end
        if abs(FR_check(i,j)) < 1e-13
            FR_check(i,j) = 0;
        end
    end
end

%% make Derivate of Fourier matrix (theta)
Dck_theta = zeros(N,N);
for j = 1:N % row
    for k = 1:N % column
        K = -(N-1)/2-1+k;
        Dck_theta(j,k) = (1i)*K/N*exp(2*(1i)*pi*(K*(j-1)/N)); % derivate of Ick
%        Dck(j,k) = 2*pi/N*(1i)*K/N*exp(2*(1i)*pi*(K*(j-1)/N)); % derivate
%        of Ick これはだめだった！
    end
end
DDFT_theta = kron(Dck_theta, eye(d));

%% make DFR_theta
DFR_theta = DDFT_theta*DFT;
for i = 1:d*N
    for j = 1:d*N
        if abs(imag(DFR_theta(i,j))) < 1e-13
            DFR_theta(i,j) = real(DFR_theta(i,j));
        end
        if abs(FR(i,j)) < 1e-13
            DFR_theta(i,j) = 0;
        end
    end
end

%% make Derivate of Fourier matrix (rho)
Dck_rho = zeros(N,N);
for j = 1:N % row
    for k = 1:N % column
        K = -(N-1)/2-1+k;
        Dck_rho(j,k) = -(1i)*K/N*exp(2*(1i)*pi*(K*(j-1)/N))*exp(-1i*K*rho); % derivate of Ick
    end
end
DDFT_rho = kron(Dck_rho, eye(d));

%% make DFR_rho
DFR_rho = DDFT_rho*DFT;
for i = 1:d*N
    for j = 1:d*N
        if abs(imag(DFR_rho(i,j))) < 1e-13
            DFR_rho(i,j) = real(DFR_rho(i,j));
        end
        if abs(FR(i,j)) < 1e-13
            DFR_rho(i,j) = 0;
        end
    end
end

end