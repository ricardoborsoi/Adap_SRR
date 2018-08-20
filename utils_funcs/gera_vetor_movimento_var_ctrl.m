function [Coord] = gera_vetor_movimento_var_ctrl(n_mov, passo, desvio)
%---- ajuda para gera_vetor_movimento_var_ctrl.m --------------------------
%
% function [Coord] = gera_vetor_movimento(n_mov, passo, desvio)
%
% Gera vetor com coordenadas de movimento global (translação). 
% Tipo de movimento: random walk com passo "passo".
%
% Parâmetros de entrada:    n_mov   = número de matrizes de movimento 
%                                     (no <= n <= nf)
%                           passo   = passo do random walk
%                           desvio  = desvio (x passos) máximo do centro 
%                                     do movimento 
%                                     (zero = não controla desvio)
%
% Parâmetros de saída:      Coord   = coordenadas (x,y) do movimento
%
% Guilherme Holsbach Costa 
% 28/06/2004
%--------------------------------------------------------------------------

if(mod(passo,1) == 0)
    if(desvio == 0)
        Coord = round((2*passo + 1) * rand(n_mov,2) - (2*passo + 1)/2);
        Coord(1,:) = [0 0];
    else
        Coord = zeros(n_mov,2);
        soma_x = 0;
        soma_y = 0;
        for i = 2:n_mov,
            flag = 1;
            while (flag),
                passo_x = round((2*passo + 1) * rand - (2*passo + 1)/2);
                soma_temp = soma_x + passo_x;
                if(abs(soma_temp) > abs(desvio*passo))
                    flag = 1;
                else
                    flag = 0;
                end
            end
            soma_x = soma_temp;
                        
            flag = 1;
            while (flag),
                passo_y = round((2*passo + 1) * rand - (2*passo + 1)/2);
                soma_temp = soma_y + passo_y;
                if(abs(soma_temp) > abs(desvio*passo))
                    flag = 1;
                else
                    flag = 0;
                end
            end
            soma_y = soma_temp;
            
            Coord(i,1) = passo_x;
            Coord(i,2) = passo_y;
        end
    end
else
    if(desvio == 0)
        Coord = 2*passo * rand(n_mov,2) - passo;
        Coord(1,:) = [0 0];
    else
        Coord = zeros(n_mov,2);
        soma_x = 0;
        soma_y = 0;
        for i = 2:n_mov,
            flag = 1;
            while (flag),
                passo_x = 2*passo * rand - passo;
                soma_temp = soma_x + passo_x;
                if(abs(soma_temp) > abs(desvio*passo))
                    flag = 1;
                else
                    flag = 0;
                end
            end
            soma_x = soma_temp;
                        
            flag = 1;
            while (flag),
                passo_y = 2*passo * rand - passo;
                soma_temp = soma_y + passo_y;
                if(abs(soma_temp) > abs(desvio*passo))
                    flag = 1;
                else
                    flag = 0;
                end
            end
            soma_y = soma_temp;
            
            Coord(i,1) = soma_x;
            Coord(i,2) = soma_y;
        end
    end
end

    

