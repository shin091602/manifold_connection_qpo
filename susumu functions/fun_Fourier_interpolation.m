function FR = fun_Fourier_interpolation(rho, p)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Dictionary
d = p("d");
N = p("N");
num_iter = p("num_iter");
num_iter = num_iter/3;

%% make DFT matrix
ck = zeros(num_iter,num_iter);
for k = 1:num_iter
    K = -(num_iter-1)/2-1+k;
    for j = 1:num_iter
        ck(k,j) = exp(-2*pi*(1i)*K*(j-1)/num_iter);
    end
end
DFT = kron(ck, eye(d)); % Discrete Fourier Transform matrix

%% make IDFT matrix
Ick_Rho = zeros(num_iter,num_iter);
for j = 1:num_iter
    for k = 1:num_iter
        K = -(num_iter-1)/2-1+k;
        Ick_Rho(j,k) = 1/exp(2*pi*(1i)*K*(j-1)/num_iter)*exp(-1i*K*rho);
    end
end
IDFT_rho = kron(Ick_Rho, eye(d)); % Inverse Discrete Fourier Transform matrix

Ick = zeros(num_iter,num_iter);
for j = 1:num_iter
    for k = 1:num_iter
        K = -(num_iter-1)/2-1+k;
        Ick(j,k) = 1/N*exp(2*pi*(1i)*K*(j-1)/num_iter);
    end
end
check_DFT = Ick*ck;
for i = 1:num_iter
    for j = 1:num_iter
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
for i = 1:d*num_iter
    for j = 1:d*num_iter
        if abs(imag(FR(i,j))) < 1e-13
            FR(i,j) = real(FR(i,j));
        end
        if abs(FR(i,j)) < 1e-13
            FR(i,j) = 0;
        end
    end
end

%% make FR_check
IDFT = kron(Ick, eye(d));
FR_check = IDFT*DFT;
for i = 1:d*num_iter
    for j = 1:d*num_iter
        if abs(imag(FR_check(i,j))) < 1e-13
            FR_check(i,j) = real(FR_check(i,j));
        end
        if abs(FR_check(i,j)) < 1e-13
            FR_check(i,j) = 0;
        end
    end
end

end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%