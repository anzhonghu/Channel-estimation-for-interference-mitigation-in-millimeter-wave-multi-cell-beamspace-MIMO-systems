%esprit approach
%esprit approach
Xa0 = sign(randn(K, N));
Xa1 = sign(randn(K, N));
Xa2 = sign(randn(K, N));
for k = 1 : K
    for l = 1 : N
        if 0 == Xa0(k,l)
            Xa0(k,l) = 1;
        else
        end
        if 0 == Xa1(k,l)
            Xa1(k,l) = 1;
        else
        end
        if 0 == Xa2(k,l)
            Xa2(k,l) = 1;
        else
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%output: U3b, spatailf3_ne, in order
beamselectforesprit_added;%KSII paper choose the beams
spasttst = zeros(3*K, 2);
spasttstforgrid = zeros(3*K, 2);
spagridsttst = zeros(K, 2);
spatailf3_ne = zeros(K, 2);
U3bb = zeros(n_arr(3, 1), 2*K);
C3_ind_sto = zeros(K, 1);
allowcount = 0;
for kx = 1 : K
    x2y2st = zeros(3, 1);
    for kx_ano = 1 : 3
        k = kx_ano + (kx-1)*3;
        spa_stt = zeros(1, 2);
        if abs(C3(k, 1)/naz(3, 1) - floor(C3(k, 1)/naz(3, 1)))<1e-2
            x_t = naz(3, 1);
            y_t = floor(C3(k, 1)/naz(3, 1));
        else
            y_t = floor(C3(k, 1)/naz(3, 1)) + 1;
            x_t = C3(k, 1) - (y_t-1) * naz(3, 1);
        end
        spa_stt(1, 1) = x_t/naz(3, 1)*2-1;
        spa_stt(1, 2) = y_t/nel(3, 1)*2-1;
        sphiaz = spa_stt(1, 1) / sqrt(1 - (spa_stt(1, 2))^2);
        tphiaz = sphiaz / sqrt(1 - sphiaz^2);
        aforypie = (spa_stt(1, 2))^2 / (cos(beta_m))^2 + tphiaz^2 * (spa_stt(1, 2))^2 - (tan(beta_m))^2;
        bforypie = 2 * height * sin(beta_m) * (1 - (spa_stt(1, 2))^2) / (cos(beta_m))^2;
        cforypie = (spa_stt(1, 2))^2 * (tan(beta_m))^2 + (spa_stt(1, 2))^2 - (tan(beta_m))^2 * (sin(beta_m))^2 - (cos(beta_m))^2 - 2 * (sin(beta_m))^2;
        cforypie = height^2 * cforypie;
        ypie1 = (-bforypie + sqrt(bforypie^2 - 4 * aforypie * cforypie)) / aforypie * 0.5;
        ypie2 = (-bforypie - sqrt(bforypie^2 - 4 * aforypie * cforypie)) / aforypie * 0.5;
        x_temp1 = (ypie1 - height * sin(beta_m)) / cos(beta_m);
        x_temp2 = (ypie2 - height * sin(beta_m)) / cos(beta_m);
        indic1 = (x_temp1 * sin(beta_m) - height * cos(beta_m)) * spa_stt(1, 2);
        indic2 = (x_temp2 * sin(beta_m) - height * cos(beta_m)) * spa_stt(1, 2);
        if ypie1 < 0  || indic1 < 0
            if ypie2 < 0 || indic2 < 0
                flagforesprit = 1;
            else
                ypie = ypie2;
            end
        else
            if ypie2 < 0 || indic2 < 0
                ypie = ypie1;
            else
                flagforesprit = 1;
            end
        end
        if 1 == flagforesprit
            flagforesprit = 0;
            x2y2 = 1e4 * Rmax^2;
        else
            x2y2 = (abs(ypie/cos(beta_m)) - height * tan(beta_m) - height * cot(beta_m))^2;
            x2y2 =  (sin(beta_m) )^2 / (spa_stt(1, 2))^2 * x2y2;
            x2y2 = x2y2 - height^2;
        end
        x2y2st(kx_ano, 1) = x2y2;
        spasttst(k, 1) = spa_stt(1, 1);
        spasttst(k, 2) = spa_stt(1, 2);
        spasttstforgrid(k, 1) = x_t;
        spasttstforgrid(k, 2) = y_t;
    end
    [~, x2y2stind] = min(x2y2st(:, 1));
    spatailf3_ne(kx, :) = spasttst((kx-1)*3 + x2y2stind, :);
    spagridsttst(kx, :) = spasttstforgrid((kx-1)*3 + x2y2stind, :);
    U3bb(:, (kx-1)*beam_numberforMS+1) = U3(:, C3((kx-1)*3 + x2y2stind, 1));
    U3bb(:, (kx-1)*beam_numberforMS+2) = U3(:, C3(3*K+((kx-1)*3 + x2y2stind-1)*(beam_numberforMS-1)+1, 1));
    C3_ind_sto(kx, 1) = (kx-1)*3 + x2y2stind;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%55
y3 = (sqrt(rho) * (H3qq * Xa0 + Hac13qq * Xa1 + Hac23qq * Xa2) + randn(n_arr(3, 1), N));
J1y3 = zeros(nq_arr(3,1), N);
U3bq1 = zeros(nq_arr(3, 1), 3*Mb);
for ny = 0 : nel(3, 1)-2
    for nx = 0 : naz(3, 1)-2
        n = ny * naz(3, 1) + 1 + nx;
        n1 = n - ny;
        J1y3(n1, :) = y3(n, :);
        U3bq1(n1, :) = U3b(n, :);
    end
end
J2y3 = zeros(nq_arr(3,1), N);
%U3bq2 = zeros(nq_arr(3, 1), Mb);
for ny = 0 : nel(3, 1)-2
    for nx = 1 : naz(3, 1)-1
        n = ny * naz(3, 1) + 1 + nx;
        n1 = n - ny - 1;
        J2y3(n1, :) = y3(n, :);
        %U3bq2(n1, :) = U3b(n, :);
    end
end
J3y3 = zeros(nq_arr(3,1), N);
%U3bq3 = zeros(nq_arr(3, 1), Mb);
for ny = 1 : nel(3, 1)-1
    for nx = 0 : naz(3, 1)-2
        n = ny * naz(3, 1) + 1 + nx;
        n1 = n - ny + 1 - naz(3, 1);
        J3y3(n1, :) = y3(n, :);
        %U3bq3(n1, :) = U3b(n, :);
    end
end
%%%%%%%%%%%%%%
spaf1 = zeros(K, 3);
spaf2 = zeros(K, 3);
spaf3 = zeros(K, 3);
spaf1_temp = zeros(3*K, 3);
spaf2_temp = zeros(3*K, 3);
spaf3_temp = zeros(3*K, 3);
spaf1_t = zeros(3*K, 3);
spaf2_t = zeros(3*K, 3);
spaf3_t = zeros(3*K, 3);
spagridsttst_t = zeros(3*K, 2);
h3 = zeros(n_arr(3, 1), 1);
Rd31 = (U3bq1' * J1y3);
Rd32 = (U3bq1' * J2y3);
Rd33 = (U3bq1' * J3y3);
R3j2 = [Rd31; Rd32] * [Rd31; Rd32]';
R3j3 = [Rd31; Rd33] * [Rd31; Rd33]';
[V12, D12] = eig(R3j2);
[~, Ind12] = sort(sum(D12), 'descend');
ExEy = [V12(1:3*K*beam_numberforMS, Ind12(1,1:3*K)),V12(3*K*beam_numberforMS+1:2*3*K*beam_numberforMS, Ind12(1,1:3*K))];
[Vexexy, Dexey] = eig(ExEy' * ExEy);
[~, Indxy] = sort(sum(Dexey));
[~, D_temp] = eig(-Vexexy(1:3*K, Indxy(1,1:3*K)) / Vexexy(3*K+1:2*3*K, Indxy(1,1:3*K)));
for k = 1 : 3*K
    spaf3_t(k, 1) = real(log(D_temp(k, k)) *1i / (2 * pi * miu));
    spagridsttst_t(k, 1) = (spaf3_t(k, 1)+1) * naz(3, 1) * 0.5;
end
spaf3_t(:,1) = sort(spaf3_t(:,1));
spagridsttst_t(:,1) = sort(spagridsttst_t(:,1));
%%%%%%%%%%%%
for k = 1 : K
    spa_com = zeros(3*K, 1);
    rp3 = rp3_store(:, k);
    nx = spagridsttst(k, 1);
    ny = spagridsttst(k, 2);
    n = (ny - 1) * naz(3, 1) + nx;
    %based on the fact that the spatail grid of each UT is not on the
    %boundary
    for kk = 1 : 3*K
        if n>2 && n<n_arr(3, 1)-1
            if (spagridsttst_t(kk,1) - spagridsttst(k,1)) < 0 && abs(spagridsttst_t(kk,1) - spagridsttst(k,1)) <= 1
                xxxx1 = sin(pi*(spatailf3_ne(k,1)-spaf3_t(kk, 1))*naz(3, 1)*0.5) / sin(pi*(spatailf3_ne(k,1)-spaf3_t(kk, 1))*0.5)...
                    /sin(pi*(spatailf3_ne(k,1)-2/naz(3, 1)-spaf3_t(kk, 1))*naz(3, 1)*0.5) * sin(pi*(spatailf3_ne(k,1)-2/naz(3, 1)-spaf3_t(kk, 1))*0.5);
                xxxx2 = abs(rp3(n, 1)) / abs(rp3(n-1, 1));
                spa_com(kk, 1) = (xxxx1-xxxx2)^2;
            else
                if (spagridsttst_t(kk,1) - spagridsttst(k,1)) > 0 && abs(spagridsttst_t(kk,1) - spagridsttst(k,1)) <= 1
                    xxxx1 = sin(pi*(spatailf3_ne(k,1)-spaf3_t(kk, 1))*naz(3, 1)*0.5) / sin(pi*(spatailf3_ne(k,1)-spaf3_t(kk, 1))*0.5)...
                        /sin(pi*(spatailf3_ne(k,1)+2/naz(3, 1)-spaf3_t(kk, 1))*naz(3, 1)*0.5) * sin(pi*(spatailf3_ne(k,1)+2/naz(3, 1)-spaf3_t(kk, 1))*0.5);
                    xxxx2 = abs(rp3(n, 1)) / abs(rp3(n+1, 1));
                    spa_com(kk, 1) = (xxxx1-xxxx2)^2;
                else
                    spa_com(kk, 1) = 1e5;
                end
            end
        else
            spa_com(kk, 1) = 1e5;
        end
    end
    [~, indforspa] = min(spa_com);
    spaf3(k, 1) = spaf3_t(indforspa, 1);
end
%%%%%%%%%%%%
[V13, D13] = eig(R3j3);
[~, Ind13] = sort(sum(D13), 'descend');
ExEy = [V13(1:3*K*beam_numberforMS, Ind13(1,1:3*K)),V13(3*K*beam_numberforMS+1:2*3*K*beam_numberforMS, Ind13(1,1:3*K))];
[Vexexy, Dexey] = eig(ExEy' * ExEy);
[~, Indxy] = sort(sum(Dexey));
[~, D_temp] = eig(-Vexexy(1:3*K, Indxy(1,1:3*K)) / Vexexy(3*K+1:2*3*K, Indxy(1,1:3*K)));
for k = 1 : 3*K
    spaf3_t(k, 2) = real(log(D_temp(k, k)) *1i / (2 * pi * miu));
    spagridsttst_t(k, 2) = (spaf3_t(k, 2)+1) * nel(3, 1) * 0.5;
end
spaf3_t(:,2) = sort(spaf3_t(:,2));
spagridsttst_t(:,2) = sort(spagridsttst_t(:,2));
%%%%%%%%%%
for k = 1 : K
    spa_com = zeros(3*K, 1);
    rp3 = rp3_store(:, k);
    nx = spagridsttst(k, 1);
    ny = spagridsttst(k, 2);
    n = (ny - 1) * naz(3, 1) + nx;
    for kk = 1 : 3*K
        if n>naz(3, 1)*2 && n<n_arr(3, 1)-naz(3, 1)*2
            if (spagridsttst_t(kk,2) - spagridsttst(k,2)) < 0 && abs(spagridsttst_t(kk,2) - spagridsttst(k,2)) <= 1
                xxxx1 = sin(pi*(spatailf3_ne(k,2)-spaf3_t(kk, 2))*nel(3, 1)*0.5) / sin(pi*(spatailf3_ne(k,2)-spaf3_t(kk, 2))*0.5)...
                    /sin(pi*(spatailf3_ne(k,2)-2/nel(3, 1)-spaf3_t(kk, 2))*nel(3, 1)*0.5) * sin(pi*(spatailf3_ne(k,2)-2/nel(3, 1)-spaf3_t(kk, 2))*0.5);
                xxxx2 = abs(rp3(n, 1)) / abs(rp3(n-naz(3, 1), 1));
                spa_com(kk, 1) = (xxxx1-xxxx2)^2;
            else
                if (spagridsttst_t(kk,2) - spagridsttst(k,2)) > 0 && abs(spagridsttst_t(kk,2) - spagridsttst(k,2)) <= 1
                    xxxx1 = sin(pi*(spatailf3_ne(k,2)-spaf3_t(kk, 2))*nel(3, 1)*0.5) / sin(pi*(spatailf3_ne(k,2)-spaf3_t(kk, 2))*0.5)...
                        /sin(pi*(spatailf3_ne(k,2)+2/nel(3, 1)-spaf3_t(kk, 2))*nel(3, 1)*0.5) * sin(pi*(spatailf3_ne(k,2)+2/nel(3, 1)-spaf3_t(kk, 2))*0.5);
                    xxxx2 = abs(rp3(n, 1)) / abs(rp3(n+naz(3, 1), 1));
                    spa_com(kk, 1) = (xxxx1-xxxx2)^2;
                else
                    spa_com(kk, 1) = 1e5;
                end
            end
        else
            spa_com(kk, 1) = 1e5;
        end
    end
    [~, indforspa] = min(spa_com);
    spaf3(k, 2) = spaf3_t(indforspa, 2);
end
%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
U3b = U3bb;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for k = 1 : K
    if abs(C3(C3_ind_sto(k, 1), 1)/naz(3, 1) - floor(C3(C3_ind_sto(k, 1), 1)/naz(3, 1)))<1e-2
        x_t = naz(3, 1);
        y_t = floor(C3(C3_ind_sto(k, 1), 1)/naz(3, 1));
    else
        y_t = floor(C3(C3_ind_sto(k, 1), 1)/naz(3, 1)) + 1;
        x_t = C3(C3_ind_sto(k, 1), 1) - (y_t-1) * naz(3, 1);
    end
    x_t = x_t/naz(3, 1)*2-1;
    y_t = y_t/nel(3, 1)*2-1;
    lambdae = exp(1i*0.5*pi*((-naz(3, 1)+1)*(x_t-spaf3(k, 1))+(-nel(3, 1)+1)*(y_t-spaf3(k, 2))));
    if x_t == spaf3(k, 1) && y_t == spaf3(k, 2)
        lambdamk = n_arr(3, 1);
    else
    end
    if x_t ~= spaf3(k, 1) && y_t == spaf3(k, 2)
        lambdamk = nel(3, 1) * (1-exp(1i*pi*naz(3, 1)*(x_t-spaf3(k, 1)))) / (1-exp(1i*pi*(x_t-spaf3(k, 1))));
    else
    end
    if x_t == spaf3(k, 1) && y_t ~= spaf3(k, 2)
        lambdamk = naz(3, 1) * (1-exp(1i*pi*nel(3, 1)*(y_t-spaf3(k, 2)))) / (1-exp(1i*pi*(y_t-spaf3(k, 2))));
    else
    end
    if x_t ~= spaf3(k, 1) && y_t ~= spaf3(k, 2)
        lambdamk = (1-exp(1i*pi*naz(3, 1)*(x_t-spaf3(k, 1)))) * (1-exp(1i*pi*nel(3, 1)*(y_t-spaf3(k, 2))))...
            / (1-exp(1i*pi*(x_t-spaf3(k, 1)))) / (1-exp(1i*pi*(y_t-spaf3(k, 2))));
    else
    end
    spaf3(k, 3) = sqrt(n_arr(3, 1)) * rp3(C3(C3_ind_sto(k, 1), 1), 1) / sqrt(rho) / lambdae / lambdamk;%strong interference may cause large error, and error is usually not small
    if abs(spaf3(k, 3)) > 1
        spaf3(k, 3) = spaf3(k, 3) / abs(spaf3(k, 3));
    else
    end
    for nx = 0 : naz(3, 1)-1
        for ny = 0 : nel(3, 1)-1
            n = ny * naz(3, 1) + 1 + nx;
            h3(n, 1) = exp(-1i * 2 * pi * miu * ((nx-0.5*(naz(3, 1)-1)) * spaf3(k, 1) + (ny-0.5*(nel(3, 1)-1)) * spaf3(k, 2)));
        end
    end
    H3_e(:, k) = h3 * spaf3(k, 3);
end
%%%%%%%%%%%%%%%%%%%%
%nmse
H_err = H3_e - H3qq;
for k = 1 : K
    err_t = 0;
    h_err_t = 0;
    for err_n = 1 : n_arr(3, 1)
        err_t = err_t + abs(H_err(err_n, k))^2;
        h_err_t = h_err_t + abs(H3qq(err_n, k))^2;
    end
    nmse(signalpo_n, 12) = nmse(signalpo_n, 12) + err_t / h_err_t;
end
%%%%%%%%%%%%%%%%%%%%%%%%
H3eb = U3b' * H3_e;
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
    capacity(signalpo_n, 12) = capacity(signalpo_n, 12) + log2(1 + sinr);
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
    capacity(signalpo_n, 12) = capacity(signalpo_n, 12) + log2(1 + sinr);
end
