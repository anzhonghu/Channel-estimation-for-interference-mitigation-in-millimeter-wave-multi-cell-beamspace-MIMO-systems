        %correaltion
        for nx = 1 : naz(3, 1)
            for ny = 1 : nel(3, 1)
                n = (ny - 1) * naz(3, 1) + nx;
                Cor3(nx, ny) = abs(rp3(n, 1));
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        temp_index3 = zeros(n_arr(3, 1), 3);
        temp_index3(:, 3) = ones(n_arr(3, 1), 1);
         %%%%%%%%%%%%%%%%%%%1111111
    for nx = 2 : naz(3, 1)-1
        for ny = 2 : nel(3, 1)-1
            n = (ny - 1) * naz(3, 1) + nx;
            if n>=in3_l && n<=in3_u && 0 ~= temp_index3(n, 3)
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
        if n>=in3_l && n<=in3_u && 0 ~= temp_index3(n, 3)
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
        if n>=in3_l && n<=in3_u && 0 ~= temp_index3(n, 3)
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
        if n>=in3_l && n<=in3_u && 0 ~= temp_index3(n, 3)
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
        if n>=in3_l && n<=in3_u && 0 ~= temp_index3(n, 3)
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
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if sum(temp_index3(:, 1)) < K
            [temp_value, temp_index_s] = sort(temp_index3(:, 3), 'descend');
            for n_com = 1 : K-sum(temp_index3(:, 1))
                temp_index3(temp_index_s(n_com, 1), 1) = 1;
                temp_index3(temp_index_s(n_com, 1), 2) = temp_value(n_com, 1);
            end
        else
            if sum(temp_index3(:, 1)) > K
                [temp_value_a, temp_index_sa] = sort(temp_index3(:, 2), 'descend');
                for n_com = K+1 : sum(temp_index3(:, 1))
                    temp_index3(temp_index_sa(n_com, 1), 1) = 0;
                    temp_index3(temp_index_sa(n_com, 1), 2) = 0;
                end
            else
            end
        end
        C3 = zeros(Mb, 2);
        nk = 0;
        index_label3 = zeros(n_arr(3, 1), 2);
        for nx = 1 : naz(3, 1)
            for ny = 1 : nel(3, 1)
                n = (ny - 1) * naz(3, 1) + nx;
                if 0 == temp_index3(n, 1)
                else
                    nk = nk + 1;
                    C3(nk, 1) = n;
                    C3(nk, 2) = temp_index3(n, 2);
                    index_label3(n, 1) = 1;
                end
            end
        end
        [~, C_index] = sort(C3(1:K,1));
        C3(1:K, :) = C3(C_index, :);
        for k = 1 : K
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
            C3(K+(k-1)*(beam_numberforMS-1)+1:K+k*(beam_numberforMS-1), 1) = index_sortan(1:(beam_numberforMS-1), 1);
            C3(K+(k-1)*(beam_numberforMS-1)+1:K+k*(beam_numberforMS-1), 2) = temp_index3(index_sortan(1:(beam_numberforMS-1), 1), 2);
            index_label3(index_sortan(1:(beam_numberforMS-1), 1), 1) = 1;
        end