%%%%%%%%%%%%%%%%%%%%%%%%%%%
in1_l3 = 42;
in1_u3 = 137;
in2_l3 = 220;
in2_u3 = 668;
in3_l3 = 900;
in3_u3 = 2662;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C3 = zeros(3*Mb, 2);
nk = 0;
index_label3 = zeros(n_arr(3, 1), 2);
rp3_store = zeros(n_arr(3, 1), K);
f_forcom = zeros(3*K, 2);
for kk = 1 : K
    rp3 = U3' * (sqrt(rho) * (H3qq(:, kk) + Hac13qq(:, kk) + Hac23qq(:, kk))   + (1/sqrt(K)) * randn(n_arr(3, 1), 1));%
    rp3_store(:, kk) = rp3;
    for nx = 1 : naz(3, 1)
        for ny = 1 : nel(3, 1)
            n = (ny - 1) * naz(3, 1) + nx;
            Cor3(nx, ny) = abs(rp3(n, 1));
        end
    end
    temp_index3 = zeros(n_arr(3, 1), 3);
    temp_index3(:, 3) = ones(n_arr(3, 1), 1);
    %%%%%%%%%%%%%%%%%%%1111111
    for nx = 2 : naz(3, 1)-1
        for ny = 2 : nel(3, 1)-1
            n = (ny - 1) * naz(3, 1) + nx;
            if 0 ~= temp_index3(n, 3)
                if Cor3(nx, ny)>Cor3(nx-1, ny) && Cor3(nx, ny)>Cor3(nx-1, ny-1) && Cor3(nx, ny)>Cor3(nx-1, ny+1) && Cor3(nx, ny)>Cor3(nx, ny-1)...
                        && Cor3(nx, ny)>Cor3(nx, ny+1) && Cor3(nx, ny)>Cor3(nx+1, ny) && Cor3(nx, ny)>Cor3(nx+1, ny-1) && Cor3(nx, ny)>Cor3(nx+1, ny+1)
                    temp_index3(n, 1) = 1;
                    temp_index3(n, 2) = Cor3(nx, ny);
                    temp_index3(n, 3) = 0;
                else
                    temp_index3(n, 3) = Cor3(nx, ny);
                end
            else
                temp_index3(n, 3) = 0;
            end
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%22222222
    nx = 1;
    for ny = 2 : nel(3, 1)-1
        n = (ny - 1) * naz(3, 1) + nx;
        if 0 ~= temp_index3(n, 3)
            if  Cor3(nx, ny)>Cor3(nx, ny-1)...
                    && Cor3(nx, ny)>Cor3(nx, ny+1) && Cor3(nx, ny)>Cor3(nx+1, ny) && Cor3(nx, ny)>Cor3(nx+1, ny-1) && Cor3(nx, ny)>Cor3(nx+1, ny+1)
                temp_index3(n, 1) = 1;
                temp_index3(n, 2) = Cor3(nx, ny);
                temp_index3(n, 3) = 0;
            else
                temp_index3(n, 3) = Cor3(nx, ny);
            end
        else
            temp_index3(n, 3) = 0;
        end
    end
    %%%%%%%%%%%%%%%%%%%33333
    nx = naz(3, 1);
    for ny = 2 : nel(3, 1)-1
        n = (ny - 1) * naz(3, 1) + nx;
        if 0 ~= temp_index3(n, 3)
            if Cor3(nx, ny)>Cor3(nx-1, ny) && Cor3(nx, ny)>Cor3(nx-1, ny-1) && Cor3(nx, ny)>Cor3(nx-1, ny+1) && Cor3(nx, ny)>Cor3(nx, ny-1)...
                    && Cor3(nx, ny)>Cor3(nx, ny+1)
                temp_index3(n, 1) = 1;
                temp_index3(n, 2) = Cor3(nx, ny);
                temp_index3(n, 3) = 0;
            else
                temp_index3(n, 3) = Cor3(nx, ny);
            end
        else
            temp_index3(n, 3) = 0;
        end
    end
    %%%%%%%%%%%%%%%%%%%4444444
    ny = 1;
    for nx = 2 : naz(3, 1)-1
        n = (ny - 1) * naz(3, 1) + nx;
        if 0 ~= temp_index3(n, 3)
            if Cor3(nx, ny)>Cor3(nx-1, ny) && Cor3(nx, ny)>Cor3(nx-1, ny+1) ...
                    && Cor3(nx, ny)>Cor3(nx, ny+1) && Cor3(nx, ny)>Cor3(nx+1, ny) && Cor3(nx, ny)>Cor3(nx+1, ny+1)
                temp_index3(n, 1) = 1;
                temp_index3(n, 2) = Cor3(nx, ny);
                temp_index3(n, 3) = 0;
            else
                temp_index3(n, 3) = Cor3(nx, ny);
            end
        else
            temp_index3(n, 3) = 0;
        end
    end
    %%%%%%%%%%%%%%55555555555
    ny = nel(3, 1);
    for nx = 2 : naz(3, 1)-1
        n = (ny - 1) * naz(3, 1) + nx;
        if 0 ~= temp_index3(n, 3)
            if Cor3(nx, ny)>Cor3(nx-1, ny) && Cor3(nx, ny)>Cor3(nx-1, ny-1) && Cor3(nx, ny)>Cor3(nx, ny-1)...
                    && Cor3(nx, ny)>Cor3(nx+1, ny) && Cor3(nx, ny)>Cor3(nx+1, ny-1)
                temp_index3(n, 1) = 1;
                temp_index3(n, 2) = Cor3(nx, ny);
                temp_index3(n, 3) = 0;
            else
                temp_index3(n, 3) = Cor3(nx, ny);
            end
        else
            temp_index3(n, 3) = 0;
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if sum(temp_index3(:, 1)) < 3*K1
        [temp_value, temp_index_s] = sort(temp_index3(:, 3), 'descend');
        for n_com = 1 : 3*K1-sum(temp_index3(:, 1))
            temp_index3(temp_index_s(n_com, 1), 1) = 1;
            temp_index3(temp_index_s(n_com, 1), 2) = temp_value(n_com, 1);
        end
    else
        if sum(temp_index3(:, 1)) > 3*K1
            [temp_value_a, temp_index_sa] = sort(temp_index3(:, 2), 'descend');
            for n_com = 3*K1+1 : sum(temp_index3(:, 1))
                temp_index3(temp_index_sa(n_com, 1), 1) = 0;
                temp_index3(temp_index_sa(n_com, 1), 2) = 0;
            end
        else
        end
    end
    for nx = 1 : naz(3, 1)
        for ny = 1 : nel(3, 1)
            n = (ny - 1) * naz(3, 1) + nx;
            if 0 == temp_index3(n, 1)
            else
                nk = nk + 1;
                C3(nk, 1) = n;
                C3(nk, 2) = temp_index3(n, 2);
                index_label3(n, 1) = 1;
                %%%%%%%%%%%%%
                f_forcom(nk, 1) = nx;
                f_forcom(nk, 2) = ny;
            end
        end
    end
end
% [~, C_index] = sort(C3(1:3*K,1));
% C3(1:3*K, :) = C3(C_index, :);
for k = 1 : 3*K
    if abs(C3(k, 1)/naz(3, 1) - floor(C3(k, 1)/naz(3, 1)))<1e-2
        x_temp = naz(3, 1);
        y_temp = floor(C3(k, 1)/naz(3, 1));
    else
        y_temp = floor(C3(k, 1)/naz(3, 1)) + 1;
        x_temp = C3(k, 1) - (y_temp-1) * naz(3, 1);
    end
    for nx = 1 : naz(3, 1)
        for ny = 1 : nel(3, 1)
            n = (ny - 1) * naz(3, 1) + nx;
            if 0 == index_label3(n, 1)
                index_label3(n, 2) = abs(nx - x_temp) + abs(ny - y_temp);
            else
                index_label3(n, 2) = n_arr(3, 1);
            end
        end
    end
    [~, index_sortan] = sort(index_label3(:, 2));
    C3(3*K+(k-1)*(beam_numberforMS-1)+1:3*K+k*(beam_numberforMS-1), 1) = index_sortan(1:(beam_numberforMS-1), 1);
    index_label3(index_sortan(1:(beam_numberforMS-1), 1), 1) = 1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
U3b = zeros(n_arr(3, 1), 3*Mb);
for k = 1 : 3*K
    U3b(:, (k-1)*beam_numberforMS+1) = U3(:, C3(k, 1));
    for nnnan = 1 : beam_numberforMS-1
        U3b(:, (k-1)*beam_numberforMS+1+nnnan) = U3(:, C3(3*K+(k-1)*(beam_numberforMS-1)+nnnan, 1));
    end
end
