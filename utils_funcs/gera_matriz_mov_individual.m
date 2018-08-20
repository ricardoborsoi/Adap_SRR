function [G] = gera_matriz_mov_individual(dx, dy, nl, nc, flag_espelha)
%--------------------------------------------------------------------------
% Gera matriz de movimento G(t).
%
% Parâmetros de entrada:    dx, dy       = vetores com as coordenadas do 
%                                          movimento de cada pixel
%                                          l(t) = l(t-1) + dy
%                                          c(t) = c(t-1) + dx
%                           nl, nc       = dimensão dos quadros a serem 
%                                          gerados
%                           flag_espelha = 0 -> considera a regiao fora da 
%                                               imagem como contendo zeros
%                                          1 -> espelha imagem ao redor 
%                                               das bordas
%                                          2 -> espelha imagem de forma
%                                               periódica
%
% Parâmetros de saída:      G        = matriz de movimento 
%
% Obs: Função não testada!
%
% Guilherme Holsbach Costa 
% 01/12/2004
%--------------------------------------------------------------------------
dx = -dx; dy = -dy;
G = sparse(nl*nc, nl*nc);    % Inicializa matriz de deslocamento

for ll = 1:nl,
    for cl = 1:nc,
        i = cl + (ll-1)*nc;
        for l = 1:nl,
            for c = 1:nc,
                j = c + (l-1)*nc;
                tempx = abs(cl - c - dx(j));
                tempy = abs(ll - l - dy(j));
                if(tempx <= 1)
                    tempx = 1 - tempx;
                    if(tempy <= 1)
                        tempy = 1 - tempy;
                        G(i,j) = tempx * tempy;
                    end
                end
                
                % imagem periódica...
                if(flag_espelha == 2)
                    tempx = abs(cl + nc - c - dx(j));
                    tempy = abs(ll + nl - l - dy(j));
                    if(tempx <= 1)
                        tempx = 1 - tempx;
                        if(tempy <= 1)
                            tempy = 1 - tempy;
                            G(i,j) = G(i,j) + tempx * tempy;
                        end
                    end
                    tempx = abs(cl - nl - c - dx(j));
                    tempy = abs(ll - nc - l - dy(j));
                    if(tempx <= 1)
                        tempx = 1 - tempx;
                        if(tempy <= 1)
                            tempy = 1 - tempy;
                            G(i,j) = G(i,j) + tempx * tempy;
                        end
                    end                    
                end
                
                % imagem espelhada...
                if(flag_espelha == 1)
                    tempx = abs(cl + nc - nc - 1 + c - dx(j));
                    tempy = abs(ll + nl - nl - 1 + l - dy(j));
                    if(tempx <= 1)
                        tempx = 1 - tempx;
                        if(tempy <= 1)
                            tempy = 1 - tempy;
                            G(i,j) = G(i,j) + tempx * tempy;
                        end
                    end
                    tempx = abs(cl - nc - nc - 1 +  c - dx(j));
                    tempy = abs(ll - nl - nl - 1 +  l - dy(j));
                    if(tempx <= 1)
                        tempx = 1 - tempx;
                        if(tempy <= 1)
                            tempy = 1 - tempy;
                            G(i,j) = G(i,j) + tempx * tempy;
                        end
                    end                    
                end
            end
        end
    end
end

%------------------------------------------
%                           flag_full= sinaliza se os lugares referentes 
%                                      às inovações devem permanecer com 
%                                      (=1) o valor da imagem anterior ou
%                                      (=0) devem conter zeros
% Funciona para movimentos inteiros
% 01/12/2004
%------------------------------------------
% for l = 1:nl,
%     for c = 1:nc,
%         i = c + (l-1)*nc;
%         l_dest = l + dy(i);
%         c_dest = c + dx(i);
%         if((l_dest <= 0) | (c_dest <= 0) | (l_dest > nl) | (c_dest > nc))
%             if(flag_full == 1)
%                 indice = i;
%                 G(i, i) = 1;
%             end
%         else
%             G(i, c_dest + (l_dest-1)*nc) = 1;
%         end
%     end
% end

