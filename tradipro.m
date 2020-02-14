H_err = Hac13qq + Hac23qq + (1/sqrt(rho)) * randn(n_arr(3, 1), K);
H3_e =  U3' * (H3qq + H_err);
D3 = zeros(Mb, 1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
diff = 1e4 * ones(n_arr(3, 1), 1);
[~, uhe3in] = sort(abs(H3_e), 'descend');
for k = 1 : K
    counttemp = 0;
    for n = 1 : n_arr(3, 1)
        if diff(uhe3in(n,k), 1) > 0
            diff(uhe3in(n,k), 1) = 0;
            counttemp = counttemp + 1;
            D3((k-1)*beam_numberforMS+counttemp, 1) = uhe3in(n,k);
        else
        end
        if counttemp < beam_numberforMS
        else
            break;
        end
    end
end
H3eb = H3_e(D3, :);
U3b = U3(:, D3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%nmse
for k = 1 : K
    err_t = 0;
    h_err_t = 0;
    for err_n = 1 : n_arr(3, 1)
        err_t = err_t + abs(H_err(err_n, k))^2;
        h_err_t = h_err_t + abs(H3qq(err_n, k))^2;
    end
    nmse(signalpo_n, 6) = nmse(signalpo_n, 6) + err_t / h_err_t;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%上行
G = H3eb / (H3eb' * H3eb + (1/rho)*eye(K));
GI = G' * U3b' * H3qq;
Xn = G' * U3b';
U3_tempfor = U3b;
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
    capacity(signalpo_n, 6) = capacity(signalpo_n, 6) + log2(1 + sinr);
end
%下行
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
    capacity(signalpo_n, 6) = capacity(signalpo_n, 6) + log2(1 + sinr);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%pilot+TDMA
%上行
G = H3eb / (H3eb' * H3eb + (1/rho)*eye(K));
GI = G' * U3b' * H3qq;
Xn = G' * U3b';
U3_tempfor = U3b;
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
    capacity(signalpo_n, 4) = capacity(signalpo_n, 4) + log2(1 + sinr);
end
%下行
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
    capacity(signalpo_n, 4) = capacity(signalpo_n, 4) + log2(1 + sinr);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% upper bound
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
diff = 1e4 * ones(n_arr(3, 1), 1);
H3_e =  U3' * (H3qq);
[~, uhe3in] = sort(abs(H3_e), 'descend');
for k = 1 : K
    counttemp = 0;
    for n = 1 : n_arr(3, 1)
        if diff(uhe3in(n,k), 1) > 0
            diff(uhe3in(n,k), 1) = 0;
            counttemp = counttemp + 1;
            D3((k-1)*beam_numberforMS+counttemp, 1) = uhe3in(n,k);
        else
        end
        if counttemp < beam_numberforMS
        else
            break;
        end
    end
end
H3eb = H3_e(D3, :);
U3b = U3(:, D3);
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
    capacity(signalpo_n, 9) = capacity(signalpo_n, 9) + log2(1 + sinr);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rhobs = rho * K;
FG = H3eb / (H3eb' * H3eb + (K/rhobs)*eye(K));
FGI = H3' * U3b * FG;
FGII1 = Hac31' * U3b * FG;%下行，假设其他两个小区内信道和参考小区内信道相同；上行，假设不同小区信道不同
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
    capacity(signalpo_n, 9) = capacity(signalpo_n, 9) + log2(1 + sinr);
end
