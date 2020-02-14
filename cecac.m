for k = 1 : K
    n_ar_in = 3;
    h = zeros(n_arr(n_ar_in, 1), 1);
    for nx = 0 : naz(n_ar_in, 1)-1
        for ny = 0 : nel(n_ar_in, 1)-1
            n = ny * naz(n_ar_in, 1) + 1 + nx;
            h(n, 1) = exp(-1i * 2 * pi * miu * ((nx-0.5*(naz(n_ar_in, 1)-1)) * spatailf3(k, 1) + (ny-0.5*(nel(n_ar_in, 1)-1)) * spatailf3(k, 2)));
        end
    end
    H3_e(:, k) = h * spatailf3(k, 3) * (beta(f_ind3(k, 4), 1)/abs(beta(f_ind3(k, 4), 1)));
end
H3eb = U3b' * H3_e;
%%%%%%%%%%%%%%%%%%%%%%%5
%nmse
H_err = H3_e - H3qq;
for k = 1 : K
    err_t = 0;
    h_err_t = 0;
    for err_n = 1 : n_arr(3, 1)
        err_t = err_t + abs(H_err(err_n, k))^2;
        h_err_t = h_err_t + abs(H3qq(err_n, k))^2;
    end
    nmse(signalpo_n, 3) = nmse(signalpo_n, 3) + err_t / h_err_t;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
G = H3eb / (H3eb' * H3eb + (1/rho)*eye(K));
GI = G' * U3b' * H3qq;
Xn = G' * U3b';
GII1 = G' * U3b' * Hac13qq;
GII2 = G' * U3b' * Hac23qq;
for k = 1 : K
    interf = 0;
    for kk = 1 : K
        if kk == k
        else
            interf = interf + abs(GI(k, kk))^2;
        end
    end
    interf = interf + GII1(k, :) * (GII1(k, :))' + GII2(k, :) * (GII2(k, :))';
    interf = rho * interf + (norm(Xn(k,:)))^2;
    signalpow = rho * abs(GI(k, k))^2;
    sinr = rho * abs(GI(k, k))^2  / interf;
    capacity(signalpo_n, 3) = capacity(signalpo_n, 3) + log2(1 + sinr);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rhobs = rho * K;
FG = H3eb / (H3eb' * H3eb + (K/rhobs)*eye(K));
FGI = H3' * U3b * FG;
FGII1 = Hac31' * U3b * FG;%假设其他两个小区内信道和参考小区内信道相同
FGII2 = Hac32' * U3b * FG;
lambdaj = zeros(K, K);
for k = 1 : K
    lambdaj(k, k) = 1 / beta(k, 1);
end
FHF = (U3b * FG)' * U3b * FG;
rhocaltemp = 0;
for k = 1 : K
    rhocaltemp = rhocaltemp + FHF(k, k) * lambdaj(k, k)^2;
end
rhocoe = sqrt(rhobs / rhocaltemp);
lambdaj = lambdaj * rhocoe;
for k = 1 : K
    interf = 0;
    for kk = 1 : K
        if kk == k
        else
            interf = interf + abs(FGI(k, kk))^2 * lambdaj(kk,kk)^2;
        end
    end
    interf = interf + FGII1(k, :) * lambdaj^2 * (FGII1(k, :))' + FGII2(k, :) * lambdaj^2 * (FGII2(k, :))';
    interf = interf + 1;
    sinr = abs(FGI(k, k))^2 * abs(lambdaj(k,k))^2  / abs(interf);
    capacity(signalpo_n, 3) = capacity(signalpo_n, 3) + log2(1 + sinr);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%5
%%%%%beamspace+TDMA
G = H3eb / (H3eb' * H3eb + (1/rho)*eye(K));
GI = G' * U3b' * H3qq;
Xn = G' * U3b';
for k = 1 : K
    interf = 0;
    for kk = 1 : K
        if kk == k
        else
            interf = interf + abs(GI(k, kk))^2;
        end
    end
    interf = rho * interf + (norm(Xn(k,:)))^2;
    signalpow = rho * abs(GI(k, k))^2;
    sinr = rho * abs(GI(k, k))^2  / interf;
    capacity(signalpo_n, 1) = capacity(signalpo_n, 1) + log2(1 + sinr);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rhobs = rho * K;
FG = H3eb / (H3eb' * H3eb + (K/rhobs)*eye(K));
FGI = H3' * U3b * FG;
lambdaj = zeros(K, K);
for k = 1 : K
    lambdaj(k, k) = 1 / beta(k, 1);
end
FHF = (U3b * FG)' * U3b * FG;
rhocaltemp = 0;
for k = 1 : K
    rhocaltemp = rhocaltemp + FHF(k, k) * lambdaj(k, k)^2;
end
rhocoe = sqrt(rhobs / rhocaltemp);
lambdaj = lambdaj * rhocoe;
for k = 1 : K
    interf = 0;
    for kk = 1 : K
        if kk == k
        else
            interf = interf + abs(FGI(k, kk))^2 * lambdaj(kk,kk)^2;
        end
    end
    interf = interf + 1;
    sinr = abs(FGI(k, k))^2 * abs(lambdaj(k,k))^2  / abs(interf);
    capacity(signalpo_n, 1) = capacity(signalpo_n, 1) + log2(1 + sinr);
end
