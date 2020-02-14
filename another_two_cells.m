%another_two_cells
base = [sqrt(3) * Rmax, Rmax; sqrt(3) * Rmax, -Rmax];
pos_ac = zeros(2*K, 3);
theta_ac = zeros(2*K, 1);
phi_ac = zeros(2*K, 1);
beta_ac = zeros(2*K, 2);
spatail_frequ_ac = zeros(2*K, 2);

pos_ac_self = zeros(2*K, 3);
theta_ac_self = zeros(2*K, 1);
phi_ac_self = zeros(2*K, 1);
beta_ac_self = zeros(2*K, 2);

Hac13 = zeros(n_arr(3, 1), K);
Hac23 = zeros(n_arr(3, 1), K);
Hac13qq = zeros(n_arr(3, 1), K);
Hac23qq = zeros(n_arr(3, 1), K);
D0f1 = zeros(K, 6);
for l = 1 : 2
    for k = 1 : K
%         pos_temp = zeros(1, 2);
        pos_temp = pos(k, 1:2);%为了简化下行计算，假设三个小区用户分布相同
%         switch l
%             case 1
%                 while norm(pos_temp) < Rmin || norm(pos_temp) > Rmax || (acos(pos_temp(1, 1) / norm(pos_temp))) < pi / 3
%                     pos_temp(1, 1) = (rand(1, 1) * 2 - 1) * Rmax;
%                     pos_temp(1, 2) = -rand(1, 1) * Rmax;
%                 end
%             case 2
%                 while norm(pos_temp) < Rmin || norm(pos_temp) > Rmax || (acos(pos_temp(1, 1) / norm(pos_temp))) < pi / 3
%                     pos_temp(1, 1) = (rand(1, 1) * 2 - 1) * Rmax;
%                     pos_temp(1, 2) = rand(1, 1) * Rmax;
%                 end
%             otherwise
%         end
        pos_ac((l-1)*K+k, 1:2) = pos_temp + base(l, :);
        pos_ac((l-1)*K+k, 3) = norm(pos_ac((l-1)*K+k, 1:2));
        pos_ac((l-1)*K+k, 3) = norm([pos_ac((l-1)*K+k, 3), height]);%distance
        phi_ac((l-1)*K+k, 1) = asin(pos_ac((l-1)*K+k, 2) / sqrt(pos_ac((l-1)*K+k, 2)^2 + (pos_ac((l-1)*K+k, 1) * cos(beta_m) + height * sin(beta_m))^2));%az
        theta_ac((l-1)*K+k, 1) = asin((pos_ac((l-1)*K+k, 1) * sin(beta_m) - height * cos(beta_m)) / pos_ac((l-1)*K+k, 3));%el
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        pos_ac_self((l-1)*K+k, 1:2) = pos_temp;
        pos_ac_self((l-1)*K+k, 3) = norm(pos_ac_self((l-1)*K+k, 1:2));
        pos_ac_self((l-1)*K+k, 3) = norm([pos_ac_self((l-1)*K+k, 3), height]);%distance
        phi_ac_self((l-1)*K+k, 1) = asin(pos_ac_self((l-1)*K+k, 2) / sqrt(pos_ac_self((l-1)*K+k, 2)^2 + (pos_ac_self((l-1)*K+k, 1) * cos(beta_m) + height * sin(beta_m))^2));%az
        theta_ac_self((l-1)*K+k, 1) = asin((pos_ac_self((l-1)*K+k, 1) * sin(beta_m) - height * cos(beta_m)) / pos_ac_self((l-1)*K+k, 3));%el
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        spatail_frequ_ac((l-1)*K+k, 1) = cos(theta_ac((l-1)*K+k, 1)) * sin(phi_ac((l-1)*K+k, 1));
        spatail_frequ_ac((l-1)*K+k, 2) = sin(theta_ac((l-1)*K+k, 1));
        f_ind_forcom_t((k-1)*3+1+l, 1) = (sin(phi_ac((l-1)*K+k, 1)) * cos(theta_ac((l-1)*K+k, 1)) + 1) * 0.5 * naz(n_ar_in, 1);
        f_ind_forcom_t((k-1)*3+1+l, 2) = (sin(theta_ac((l-1)*K+k, 1)) + 1) * 0.5 * nel(n_ar_in, 1);
    end
end
index_fortwocells = zeros(K, 2);
index_for_twocells_selection = zeros(K, 2);
for k = 1 : K
    distancefortwocells = zeros(K*K, 1);
    for kk = 1 : K
        for kkk = 1 : K
            disfortwo1 = (f_ind_forcom((k-1)*3+1, 1) - f_ind_forcom_t((kk-1)*3+1+1, 1))^2 + (f_ind_forcom((k-1)*3+1, 2) - f_ind_forcom_t((kk-1)*3+1+1, 2))^2;
            disfortwo2 = (f_ind_forcom((k-1)*3+1, 1) - f_ind_forcom_t((kkk-1)*3+1+2, 1))^2 + (f_ind_forcom((k-1)*3+1, 2) - f_ind_forcom_t((kkk-1)*3+1+2, 2))^2;
            disfortwo3 = (f_ind_forcom_t((kkk-1)*3+1+2, 1) - f_ind_forcom_t((kk-1)*3+1+1, 1))^2 + (f_ind_forcom_t((kkk-1)*3+1+2, 2) - f_ind_forcom_t((kk-1)*3+1+1, 2))^2;
            if disfortwo1 < 6 || disfortwo2 < 6 || disfortwo3 < 6
                distancefortwocells((kk-1)*K+kkk, 1) = 0;
            else
                distancefortwocells((kk-1)*K+kkk, 1) = disfortwo1 + disfortwo2 + disfortwo3;
            end
        end
    end
    [~, indextotal] = sort(distancefortwocells, 'descend');
    for twocellcount = 1 : K*K
        if rem(indextotal(twocellcount,1), K) < 1e-2
            kkk = K;
            kk = indextotal(twocellcount,1) / K;
        else
            kk = floor(indextotal(twocellcount,1) / K) + 1;
            kkk = rem(indextotal(twocellcount,1), K);
        end
        if index_for_twocells_selection(kk, 1) < 1e-2 && index_for_twocells_selection(kkk, 2) < 1e-2
            break;
        else
        end
    end
    index_for_twocells_selection(kk, 1) = 1;
    index_for_twocells_selection(kkk, 2) = 1;
    index_fortwocells(k, 1) = kk;
    index_fortwocells(k, 2) = kkk;
end
for l = 1 : 2
    for kx = 1 : K
        k = index_fortwocells(kx, l);
        f_ind_forcom((kx-1)*3+1+l, 1) = f_ind_forcom_t((k-1)*3+1+l, 1);
        f_ind_forcom((kx-1)*3+1+l, 2) = f_ind_forcom_t((k-1)*3+1+l, 2);
        Aaz = -min(12 * phi_ac((l-1)*K+k, 1)^2 / (70/180*pi)^2, 25);
        Ael = -min(12 * theta_ac((l-1)*K+k, 1)^2 / (7/180*pi)^2, 20);
        D0 = -min(-Aaz-Ael, 25);
        D0 = 10^(D0*0.1);
        D0f1(k, (l-1)*3+1) = Aaz;
        D0f1(k, (l-1)*3+2) = Ael;
        D0f1(k, (l-1)*3+3) = D0;
        beta_ac((l-1)*K+k, 1) = sqrt(D0 * lambda^2 / (16 * pi^2 * pos_ac((l-1)*K+k, 3)^2)) * exp(1i * rand(1,1) * 2 * pi);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        Aaz = -min(12 * (phi_ac_self((l-1)*K+k, 1)-pi*2/3)^2 / (70/180*pi)^2, 25);
        Ael = -min(12 * theta_ac_self((l-1)*K+k, 1)^2 / (7/180*pi)^2, 20);
        D0 = -min(-Aaz-Ael, 25);
        D0 = 10^(D0*0.1);
        beta_ac_self((l-1)*K+k, 1) = sqrt(D0 * lambda^2 / (16 * pi^2 * pos_ac_self((l-1)*K+k, 3)^2)) * exp(1i * rand(1,1) * 2 * pi);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        n_ar_in = 3;
        h = zeros(n_arr(n_ar_in, 1), 1);
        Hh = zeros(n_arr(n_ar_in, 1), 1);
        Hhqq = zeros(n_arr(n_ar_in, 1), 1);
        for path_i = 1 : path_num
            if 1 == path_i
                if (rand(1,1) - 0.1) < 0
                    beta_s(path_i, 1) = 0;
                else
                    beta_s(path_i, 1) = beta_ac((l-1)*K+k, 1);
                end
                for nx = 0 : naz(n_ar_in, 1)-1
                    for ny = 0 : nel(n_ar_in, 1)-1
                        n = ny * naz(n_ar_in, 1) + 1 + nx;
                        h(n, 1) = exp(-1i * 2 * pi * miu * ((nx-0.5*(naz(n_ar_in, 1)-1)) * sin(phi_ac((l-1)*K+k, 1)) * cos(theta_ac((l-1)*K+k, 1)) + (ny-0.5*(nel(n_ar_in, 1)-1)) * sin(theta_ac((l-1)*K+k, 1))));
                    end
                end
            else
                beta_s(path_i, 1) = beta_ac((l-1)*K+k, 1) * 10^((-rand(1,1) * 5 - 15)*0.1);
                theta_ac_rand = theta_ac((l-1)*K+k, 1)+(rand(1,1)-0.5)*6*pi/180;
                phi_ac_rand = phi_ac((l-1)*K+k, 1)+(rand(1,1)-0.5)*20*pi/180;
                for nx = 0 : naz(n_ar_in, 1)-1
                    for ny = 0 : nel(n_ar_in, 1)-1
                        n = ny * naz(n_ar_in, 1) + 1 + nx;
                        h(n, 1) = exp(-1i * 2 * pi * miu * ((nx-0.5*(naz(n_ar_in, 1)-1)) * sin(phi_ac_rand) ...
                            * cos(theta_ac_rand) + (ny-0.5*(nel(n_ar_in, 1)-1)) * sin(theta_ac_rand)));
                    end
                end
            end
            Hh = Hh + h * beta_s(path_i, 1);
        end
        Hhqq = Hh / max(abs(beta_s(:, 1)));
        %%%%%%%%%%%%%%%%%%%%%%%%%
        switch l
            case 1
                Hac13(:, kx) = Hh;
                Hac13qq(:, kx) = Hhqq;
            case 2
                Hac23(:, kx) = Hh;
                Hac23qq(:, kx) = Hhqq;
            otherwise
        end
    end
end