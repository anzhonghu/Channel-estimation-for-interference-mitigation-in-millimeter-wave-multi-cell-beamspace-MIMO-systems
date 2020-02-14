clear;
close all;
naz = [7; 16; 32];
nel = [36; 76; 151];
n_arr = naz .* nel;
nq_arr = (naz-1) .* (nel-1);
K = 50;
Mb = 2 * K;
beam_numberforMS = floor(Mb/K);
Rmin = 10;
Rmax = 100;
rho_store = [1e-2; 10^(-1); 1e0; 10^(1); 1e2;];
height = 10;
f = 80*10^9;%1G bandwidth
lambda = 3 * 10^8 / f;
miu = 0.5;
Tc = 500;
N = Tc-K;
N_1 = Tc - 1;
Spad = 5;
K1 = 1;
path_num = 4;
betam_min = atan(height/Rmax);
betam_max = atan(height/Rmin);
% beta_m = 0.5*(betam_min+betam_max);%larger than 0.5*(betam_min+betam_max)
beta_m = 10.3 * pi / 180;%-phi_m in the EL paper
Nite = 1e3;%1e3;
c_km = zeros(K, beam_numberforMS);
el_in = zeros(Mb, 1);
az_in = zeros(Mb, 1);
capacity = zeros(length(rho_store),  12);
nmse = zeros(length(rho_store),  12);
esti_err = zeros(length(rho_store),  3);
angle = zeros(1, 2);
flagforesprit = 0;
U1 = zeros(n_arr(1, 1), n_arr(1, 1));
U2 = zeros(n_arr(2, 1), n_arr(2, 1));
U3 = zeros(n_arr(3, 1), n_arr(3, 1));
for nx = 1 : naz(3, 1)
    for ny = 1 : nel(3, 1)
        angle(1, 2) = (-1+2*ny/nel(3, 1));%el
        angle(1, 1) = (-1+2*nx/naz(3, 1));%az
        n = (ny - 1) * naz(3, 1) + nx;
        for mx = 0 : naz(3, 1)-1
            for my = 0 : nel(3, 1)-1
                m = my * naz(3, 1) + 1 + mx;
                U3(m, n) = exp(-1i * 2 * pi * miu * ((mx-0.5*(naz(3, 1)-1)) * angle(1, 1) + (my-0.5*(nel(3, 1)-1)) * angle(1, 2))) / sqrt(n_arr(3, 1));
            end
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%
in1_l = 42;
in1_u = 132;
in2_l = 220;
in2_u = 648;
in3_l = 900;
in3_u = 2593;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



H3_e = zeros(n_arr(3, 1), K);

H3eb = zeros(Mb, K);
H3qq = zeros(n_arr(3, 1), K);
Cor3 = zeros(naz(3, 1), nel(3, 1));
for signalpo_n = 1 : length(rho_store)
    rho = rho_store(signalpo_n, 1);
    for ii = 1 : Nite
        pos = zeros(K, 3);
        pos31 = zeros(K, 3);
        pos32 = zeros(K, 3);
        theta = zeros(K, 1);
        theta_rand_s = zeros(path_num, 1);
        phi_rand_s = zeros(path_num, 1);
        beta_s = zeros(path_num, 1);
        theta_rand_store = zeros(path_num, K);
        phi_rand_store = zeros(path_num, K);
        beta_store = zeros(path_num, K);
        phi = zeros(K, 1);
        theta31 = zeros(K, 1);
        theta32 = zeros(K, 1);
        phi31 = zeros(K, 1);
        phi32 = zeros(K, 1);
        f_ind1 = zeros(K, 4);
        f_ind2 = zeros(K, 4);
        f_ind3 = zeros(K, 4);
        f_ind_forcom_t = zeros(3*K, 2);
        f_ind_forcom = zeros(3*K, 2);
        spatail_frequ = zeros(K, 2);
        beta = zeros(K, 2);
        beta_aa = zeros(K, 2);
        beta31 = zeros(K, 2);
        beta32 = zeros(K, 2);
        beam_in = zeros(Mb, 1);
        H3 = zeros(n_arr(3, 1), K);
        Hac31 = zeros(n_arr(3, 1), K);
        Hac32 = zeros(n_arr(3, 1), K);
        for k = 1 : K
            pos_temp = zeros(1, 2);
            while norm(pos_temp) < Rmin || norm(pos_temp) > Rmax || abs(atan(pos_temp(1, 2) / pos_temp(1, 1))) > pi / 3
                pos_temp(1, 1) = rand(1, 1) * Rmax;
                pos_temp(1, 2) = (rand(1, 1) * 2 - 1) * Rmax;
            end
            pos(k, 1:2) = pos_temp;
            pos(k, 3) = norm(pos_temp);
            pos(k, 3) = norm([pos(k, 3), height]);%distance
            phi(k, 1) = asin(pos_temp(1, 2) / sqrt(pos_temp(1, 2)^2 + (pos_temp(1, 1) * cos(beta_m) + height * sin(beta_m))^2));%az
            theta(k, 1) = asin((pos_temp(1, 1) * sin(beta_m) - height * cos(beta_m)) / pos(k, 3));%el
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            Aaz = -min(12 * phi(k, 1)^2 / (70/180*pi)^2, 25);
            Ael = -min(12 * theta(k, 1)^2 / (7/180*pi)^2, 20);
            D0 = -min(-Aaz-Ael, 25);
            D0 = 10^(D0*0.1);
            beta(k, 1) = sqrt(D0 * lambda^2 / (16 * pi^2 * pos(k, 3)^2)) * exp(1i * rand(1,1) * 2 * pi);
            for path_i = 1 : path_num
                if 1 == path_i
                    if (rand(1,1) - 0.1) < 0
                        beta_store(path_i, k) = 0;
                    else
                        beta_store(path_i, k) = beta(k, 1);
                    end
                    theta_rand_store(path_i, k) = theta(k, 1);
                    phi_rand_store(path_i, k) = phi(k, 1);
                else
                    beta_store(path_i, k) = beta(k, 1) * 10^((-rand(1,1) * 5 - 15)*0.05);
                    theta_rand_store(path_i, k) = theta(k, 1)+(rand(1,1)-0.5)*6*pi/180;
                    phi_rand_store(path_i, k) = phi(k, 1)+(rand(1,1)-0.5)*20*pi/180;
                end
            end
            [~,beta_index_a] = max(abs(beta_store(:, k)));
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            pos31(k, 1:2) = pos_temp + [sqrt(3) * Rmax, Rmax];
            pos31(k, 3) = norm(pos31(k, 1:2));
            pos31(k, 3) = norm([pos31(k, 3), height]);%distance
            phi31(k, 1) = asin(pos31(k, 2) / sqrt(pos31(k, 2)^2 + (pos31(k, 1) * cos(beta_m) + height * sin(beta_m))^2));%az
            theta31(k, 1) = asin((pos31(k, 1) * sin(beta_m) - height * cos(beta_m)) / pos31(k, 3));%el
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            pos32(k, 1:2) = pos_temp + [sqrt(3) * Rmax, -Rmax];
            pos32(k, 3) = norm(pos32(k, 1:2));
            pos32(k, 3) = norm([pos32(k, 3), height]);%distance
            phi32(k, 1) = asin(pos32(k, 2) / sqrt(pos32(k, 2)^2 + (pos32(k, 1) * cos(beta_m) + height * sin(beta_m))^2));%az
            theta32(k, 1) = asin((pos32(k, 1) * sin(beta_m) - height * cos(beta_m)) / pos32(k, 3));%el
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            spatail_frequ(k, 1) = cos(theta_rand_store(beta_index_a, k)) * sin(phi_rand_store(beta_index_a, k));
            spatail_frequ(k, 2) = sin(theta_rand_store(beta_index_a, k));
            n_ar_in = 3;
            f_ind3(k, 1) = (spatail_frequ(k, 1) + 1) * 0.5 * naz(n_ar_in, 1);
            f_ind3(k, 2) = (spatail_frequ(k, 2) + 1) * 0.5 * nel(n_ar_in, 1);
            f_ind3(k, 3) = (f_ind3(k, 2) - 1) * naz(n_ar_in, 1) + f_ind3(k, 1);
            f_ind_forcom_t((k-1)*3+1, 1) = f_ind3(k, 1);
            f_ind_forcom_t((k-1)*3+1, 2) = f_ind3(k, 2);
        end
        [~, f_ind3(:, 4)] = sort(f_ind3(:, 3));
        spatial_forcom = zeros(K, 2);
        D0forcom = zeros(K, 3);
        for kx = 1 : K
            n_ar_in =  3;
            k = f_ind3(kx, 4);
            f_ind_forcom((kx-1)*3+1, 1) = f_ind_forcom_t((k-1)*3+1, 1);
            f_ind_forcom((kx-1)*3+1, 2) = f_ind_forcom_t((k-1)*3+1, 2);
            Aaz = -min(12 * phi(k, 1)^2 / (70/180*pi)^2, 25);
            Ael = -min(12 * theta(k, 1)^2 / (7/180*pi)^2, 20);
            D0 = -min(-Aaz-Ael, 25);
            D0 = 10^(D0*0.1);
            D0forcom(kx, 1) = Aaz;
            D0forcom(kx, 2) = Ael;
            D0forcom(kx, 3) = D0;
            beta(k, 1) = sqrt(D0 * lambda^2 / (16 * pi^2 * pos(k, 3)^2)) * exp(1i * rand(1,1) * 2 * pi);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            Aaz = -min(12 * phi31(k, 1)^2 / (70/180*pi)^2, 25);
            Ael = -min(12 * theta31(k, 1)^2 / (7/180*pi)^2, 20);
            D0 = -min(-Aaz-Ael, 25);
            D0 = 10^(D0*0.1);
            beta31(k, 1) = sqrt(D0 * lambda^2 / (16 * pi^2 * pos31(k, 3)^2)) * exp(1i * rand(1,1) * 2 * pi);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            Aaz = -min(12 * phi32(k, 1)^2 / (70/180*pi)^2, 25);
            Ael = -min(12 * theta32(k, 1)^2 / (7/180*pi)^2, 20);
            D0 = -min(-Aaz-Ael, 25);
            D0 = 10^(D0*0.1);
            beta32(k, 1) = sqrt(D0 * lambda^2 / (16 * pi^2 * pos32(k, 3)^2)) * exp(1i * rand(1,1) * 2 * pi);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            h = zeros(n_arr(n_ar_in, 1), 1);
            for path_i = 1 : path_num
                if 1 == path_i
                    for nx = 0 : naz(n_ar_in, 1)-1
                        for ny = 0 : nel(n_ar_in, 1)-1
                            n = ny * naz(n_ar_in, 1) + 1 + nx;
                            h(n, 1) = exp(-1i * 2 * pi * miu * ((nx-0.5*(naz(n_ar_in, 1)-1)) * sin(phi(k, 1)) * cos(theta(k, 1)) + (ny-0.5*(nel(n_ar_in, 1)-1)) * sin(theta(k, 1))));
                        end
                    end
                else
                    for nx = 0 : naz(n_ar_in, 1)-1
                        for ny = 0 : nel(n_ar_in, 1)-1
                            n = ny * naz(n_ar_in, 1) + 1 + nx;
                            h(n, 1) = exp(-1i * 2 * pi * miu * ((nx-0.5*(naz(n_ar_in, 1)-1)) * sin(phi_rand_store(path_i, k)) * cos(theta_rand_store(path_i, k)) + (ny-0.5*(nel(n_ar_in, 1)-1)) * sin(theta_rand_store(path_i, k))));
                        end
                    end
                end
                H3(:, kx) = H3(:, kx) + h * beta_store(path_i, k);
            end
            [~,beta_index_a] = max(abs(beta_store(:, k)));
            beta_aa(k, 1) = beta_store(beta_index_a, k);
            H3qq(:, kx) = H3(:, kx) / abs(beta_aa(k, 1));
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            for path_i = 1 : path_num
                if 1 == path_i
                    if (rand(1,1) - 0.1) < 0
                        beta_s(path_i, 1) = 0;
                    else
                        beta_s(path_i, 1) = beta31(k, 1);
                    end
                    for nx = 0 : naz(n_ar_in, 1)-1
                        for ny = 0 : nel(n_ar_in, 1)-1
                            n = ny * naz(n_ar_in, 1) + 1 + nx;
                            h(n, 1) = exp(-1i * 2 * pi * miu * ((nx-0.5*(naz(n_ar_in, 1)-1)) * sin(phi31(k, 1)) * cos(theta31(k, 1)) + (ny-0.5*(nel(n_ar_in, 1)-1)) * sin(theta31(k, 1))));
                        end
                    end
                else
                    beta_s(path_i, 1) = beta31(k, 1) * 10^((-rand(1,1) * 5 - 15)*0.1);
                    theta_ac_rand = theta31(k, 1)+(rand(1,1)-0.5)*6*pi/180;
                    phi_ac_rand = phi31(k, 1)+(rand(1,1)-0.5)*20*pi/180;
                    for nx = 0 : naz(n_ar_in, 1)-1
                        for ny = 0 : nel(n_ar_in, 1)-1
                            n = ny * naz(n_ar_in, 1) + 1 + nx;
                             h(n, 1) = exp(-1i * 2 * pi * miu * ((nx-0.5*(naz(n_ar_in, 1)-1)) * sin(phi_ac_rand) * cos(theta_ac_rand) + (ny-0.5*(nel(n_ar_in, 1)-1)) * sin(theta_ac_rand)));
                        end
                    end
                end
                Hac31(:, kx) = Hac31(:, kx) + h * beta_s(path_i, 1);
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            for path_i = 1 : path_num
                if 1 == path_i
                    if (rand(1,1) - 0.1) < 0
                        beta_s(path_i, 1) = 0;
                    else
                        beta_s(path_i, 1) = beta32(k, 1);
                    end
                    for nx = 0 : naz(n_ar_in, 1)-1
                        for ny = 0 : nel(n_ar_in, 1)-1
                            n = ny * naz(n_ar_in, 1) + 1 + nx;
                            h(n, 1) = exp(-1i * 2 * pi * miu * ((nx-0.5*(naz(n_ar_in, 1)-1)) * sin(phi32(k, 1)) * cos(theta32(k, 1)) + (ny-0.5*(nel(n_ar_in, 1)-1)) * sin(theta32(k, 1))));
                        end
                    end
                else
                    beta_s(path_i, 1) = beta32(k, 1) * 10^((-rand(1,1) * 5 - 15)*0.1);
                    theta_ac_rand = theta32(k, 1)+(rand(1,1)-0.5)*6*pi/180;
                    phi_ac_rand = phi32(k, 1)+(rand(1,1)-0.5)*20*pi/180;
                    for nx = 0 : naz(n_ar_in, 1)-1
                        for ny = 0 : nel(n_ar_in, 1)-1
                            n = ny * naz(n_ar_in, 1) + 1 + nx;
                             h(n, 1) = exp(-1i * 2 * pi * miu * ((nx-0.5*(naz(n_ar_in, 1)-1)) * sin(phi_ac_rand) * cos(theta_ac_rand) + (ny-0.5*(nel(n_ar_in, 1)-1)) * sin(theta_ac_rand)));
                        end
                    end
                end
                Hac32(:, kx) = Hac32(:, kx) + h * beta_s(path_i, 1);
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            spatial_forcom(kx, 1) = spatail_frequ(k, 1);
            spatial_forcom(kx, 2) = spatail_frequ(k, 2);
        end
        another_two_cells;
        spaxx = [spatail_frequ(:, 1); spatail_frequ_ac(:, 1)];
        spayy = [spatail_frequ(:, 2); spatail_frequ_ac(:, 2)];
        [~, az_ind] = sort([spatail_frequ(:, 1); spatail_frequ_ac(:, 1)]);
        [~, el_ind] = sort([spatail_frequ(:, 2); spatail_frequ_ac(:, 2)]);
        rp3 = (U3' * (sqrt(rho) * (sum(H3qq, 2)  + sum(Hac13qq, 2) + sum(Hac23qq, 2)) + randn(n_arr(3, 1), 1)));
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %traditional processing
        tradipro;
        %%%%%%%%%%%%%%%%%%%%%%
        corbes;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %generate beams, DOA estimation
        gbdoae;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %channel estimation, capacity calculation
        cecac;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
        %proposed approach
        esprit_c3_added;
        disp([signalpo_n, ii])
    end
end
capacity = capacity / Nite;
nmse = nmse / Nite;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
capacity(:, 1) = capacity(:, 1) * N_1 * 0.5 / 3;%%%%%%%%TDMA+[22]
capacity(:, 3) = capacity(:, 3) * N_1 * 0.5;
capacity(:, 7:9) = capacity(:, 7:9) * Tc * 0.5;
capacity(:, 4) = capacity(:, 4) * N * 0.5 / 3;%%%%%%%%TDMA+[21]
capacity(:, 6) = capacity(:, 6) * N * 0.5;
capacity(:, 10:12) = capacity(:, 10:12) * N * 0.5;
capacity = capacity / Tc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h = figure;
set(h,'PaperType','A4');
axes('FontSize',16);
subplot(1,2,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot(10 * log10(rho_store), capacity(1:length(rho_store), 3), 'k--d','LineWidth',2,'MarkerSize',14)
hold on
plot(10 * log10(rho_store), capacity(1:length(rho_store), 1), 'k--^','LineWidth',2,'MarkerSize',14)
plot(10 * log10(rho_store), capacity(1:length(rho_store), 6), 'k--x','LineWidth',2,'MarkerSize',10)
plot(10 * log10(rho_store), capacity(1:length(rho_store), 4), 'k--v','LineWidth',2,'MarkerSize',10)
plot(10 * log10(rho_store), capacity(1:length(rho_store), 12), 'k-o','LineWidth',2,'MarkerSize',10)
plot(10 * log10(rho_store), capacity(1:length(rho_store), 9), 'k-s','LineWidth',2,'MarkerSize',10)
xlim([min(10 * log10(rho_store)), max(10 * log10(rho_store))])
ylim([0,500])
le = legend('Beamspace based [28]','TDMA+[28]', 'Pilot based [27]','TDMA+[27]','Proposed','Perfect', 'Location', 'northwest');
set(le,'Fontsize',16,'Fontname','Times')
set(gca,'XTick',10 * log10(rho_store))
xlabel('Transmission SNR (dB)','Fontsize',20,'Fontname','Times')
ylabel('Spectral efficiency (bps/Hz)','Fontsize',20,'Fontname','Times')
grid on

subplot(1,2,2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot(10 * log10(rho_store), nmse(1:length(rho_store), 3), 'k--d','LineWidth',2,'MarkerSize',14)
hold on
plot(10 * log10(rho_store), nmse(1:length(rho_store), 6), 'k--x','LineWidth',2,'MarkerSize',10)
plot(10 * log10(rho_store), nmse(1:length(rho_store), 12), 'k-o','LineWidth',2,'MarkerSize',10)
xlim([min(10 * log10(rho_store)), max(10 * log10(rho_store))])
% ylim([0,500])
le = legend('Beamspace based [28]', 'Pilot based [27]','Proposed', 'Location', 'northwest');
set(le,'Fontsize',16,'Fontname','Times')
set(gca,'XTick',10 * log10(rho_store))
xlabel('Transmission SNR (dB)','Fontsize',20,'Fontname','Times')
ylabel('NMSE','Fontsize',20,'Fontname','Times')
grid on
print(h,'-dpdf','fig_snr')
