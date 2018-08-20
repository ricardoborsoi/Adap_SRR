function [x] = move_imagem(x, dx, dy, nl, nc, flag_espelha)
%--------------------------------------------------------------------------
% function [x] = move_imagem(x, dx, dy, nl, nc, flag_espelha)
%
% Executa movimento dos pixels de uma imagem [x_out = G*x_in].
%
%  |1 2 3|  dx = [1, ..., 1]   |5 6 x|
%  |4 5 6|  ================>  |8 9 x| 
%  |7 8 9|  dy = [1, ..., 1]   |x x x|
%
% Em que:
% x_out = x_in + dx ; y_out = y_in + dy
%
% Parâmetros de entrada:    x        = imagem de entrada (representação
%                                      lexicográfica ou não)
%                           dx, dy   = vetores (ou matrizes) com as 
%                                      coordenadas do movimento de cada 
%                                      pixel
%                           nl, nc   = dimensão dos quadros a serem gerados
%                           flag_espelha = 0 -> zero padding 
%                                          1 -> lacunas preenchidas com o 
%                                               valor do pixel da imagem
%                                               antes do deslocamento
%
% Parâmetros de saída:      x        = imagem de saída = imagem de entrada 
%                                      deslocada de dx e dy
%
% Obs: Função não implementada 
%
% Guilherme Holsbach Costa 
% 27/07/2006
%--------------------------------------------------------------------------

flag_lexico = 0;                    % indica se a imagem de entrada foi 
                                    % passada na forma lexicográfica
if(min(size(x)) == 1)               
    x = ilexico(x, nl, nc);         %
    flag_lexico = 1;                % indica que a imagem de entrada foi 
                                    % passada na forma lexicográfica
    dx = ilexico(dx, nl, nc);
    dy = ilexico(dx, nl, nc);
end

x_ = zeros(nl,nc);
mapa = zeros(nl,nc);
for l = 1:nl,
    for c = 1:nc,
        destino_c = c + dx(l,c);
        destino_l = l + dy(l,c);
        if((destino_c > 0)&(destino_c <= nc)&(destino_l > 0)&(destino_l <= nl))
            if(mod(destino_c,1) == 0)
                if(mod(destino_l,1) == 0)
                    if((destino_l > 0)&(destino_l <= nl) & (destino_c > 0)&(destino_c <= nc))
                        x_(destino_l, destino_c) = x_(destino_l, destino_c) + x(l,c);
                        mapa(destino_l, destino_c) = mapa(destino_l, destino_c) + 1;
                    end
                else
                    temp = sign(destino_l)*mod(abs(destino_l),1);
                    if((destino_l - temp > 0)&(destino_l - temp <= nl))
                        x_(destino_l - temp, destino_c) = x_(destino_l - temp, destino_c) + (1 - mod(abs(destino_l),1)) * x(l,c);
                        mapa(destino_l - temp, destino_c) = mapa(destino_l - temp, destino_c) + (1 - mod(abs(destino_l),1));
                    end
                    if((destino_l - mod(destino_l,1) + sign(destino_l) > 0)&(destino_l - mod(destino_l,1) + sign(destino_l) <= nl))
                        x_(destino_l - mod(destino_l,1) + sign(destino_l), destino_c) = x_(destino_l - mod(destino_l,1) + sign(destino_l), destino_c) + mod(abs(destino_l),1) * x(l,c);
                        mapa(destino_l - mod(destino_l,1) + sign(destino_l), destino_c) = mapa(destino_l - mod(destino_l,1) + sign(destino_l), destino_c) + mod(abs(destino_l),1);
                    end
                end
            else
                if(mod(destino_c,1) == 0)
                    temp = sign(destino_c)*mod(abs(destino_c),1);
                    if((destino_c - temp > 0)&(destino_c - temp <= nc))
                        x_(destino_l, destino_c - temp) = x_(destino_l, destino_c - temp) + (1 - mod(abs(destino_c),1)) * x(l,c);
                        mapa(destino_l, destino_c - temp) = mapa(destino_l, destino_c - temp) + (1 - mod(abs(destino_c),1));
                    end
                    if((destino_c - mod(destino_c,1) + sign(destino_c) > 0)&(destino_c - mod(destino_c,1) + sign(destino_c) <= nl))
                        x_(destino_l, destino_c - mod(destino_c,1) + sign(destino_c)) = x_(destino_l, destino_c - mod(destino_c,1) + sign(destino_c)) + mod(abs(destino_c),1) * x(l,c);
                        mapa(destino_l, destino_c - mod(destino_c,1) + sign(destino_c)) = mapa(destino_l, destino_c - mod(destino_c,1) + sign(destino_c)) + mod(abs(destino_c),1);
                    end
                else
                    temp_c = sign(destino_c)*mod(abs(destino_c),1);
                    temp_l = sign(destino_l)*mod(abs(destino_l),1);
                    if((destino_l - temp_l > 0)&(destino_l - temp_l <= nl) & (destino_c - temp_c > 0)&(destino_c - temp_c <= nc))
                        x_(destino_l - temp_l, destino_c - temp_c) = x_(destino_l - temp_l, destino_c - temp_c) + (1 - mod(abs(destino_l),1)) * (1 - mod(abs(destino_c),1)) * x(l,c);
                        mapa(destino_l - temp_l, destino_c - temp_c) = mapa(destino_l - temp_l, destino_c - temp_c) + (1 - mod(abs(destino_l),1)) * (1 - mod(abs(destino_c),1));
                    end
                    if((destino_l - temp_l > 0)&(destino_l - temp_l <= nl) & (destino_c - mod(destino_c,1) + sign(destino_c) > 0)&(destino_c - mod(destino_c,1) + sign(destino_c) <= nc))
                        x_(destino_l - temp_l, destino_c - mod(destino_c,1) + sign(destino_c)) = x_(destino_l - temp_l, destino_c - mod(destino_c,1) + sign(destino_c)) + (1 - mod(abs(destino_l),1)) * mod(abs(destino_c),1) * x(l,c);
                        mapa(destino_l - temp_l, destino_c - mod(destino_c,1) + sign(destino_c)) = mapa(destino_l - temp_l, destino_c - mod(destino_c,1) + sign(destino_c)) + (1 - mod(abs(destino_l),1)) * mod(abs(destino_c),1);
                    end
                    if((destino_c - temp_c > 0)&(destino_c - temp_c <= nc) & (destino_l - mod(destino_l,1) + sign(destino_l) > 0)&(destino_l - mod(destino_l,1) + sign(destino_l) <= nl))
                        x_(destino_l - mod(destino_l,1) + sign(destino_l), destino_c - temp_c) = x_(destino_l - mod(destino_l,1) + sign(destino_l), destino_c - temp_c) + mod(abs(destino_l),1) * (1 - mod(abs(destino_c),1)) * x(l,c);
                        mapa(destino_l - mod(destino_l,1) + sign(destino_l), destino_c - temp_c) = mapa(destino_l - mod(destino_l,1) + sign(destino_l), destino_c - temp_c) + mod(abs(destino_l),1) * (1 - mod(abs(destino_c),1));
                    end
                    if((destino_l - mod(destino_l,1) + sign(destino_l) > 0)&(destino_l - mod(destino_l,1) + sign(destino_l) <= nl) & (destino_c - mod(destino_c,1) + sign(destino_c) > 0)&(destino_c - mod(destino_c,1) + sign(destino_c) <= nc))
                        x_(destino_l - mod(destino_l,1) + sign(destino_l), destino_c - mod(destino_c,1) + sign(destino_c)) = x_(destino_l - mod(destino_l,1) + sign(destino_l), destino_c - mod(destino_c,1) + sign(destino_c)) + mod(abs(destino_l),1) * mod(abs(destino_c),1) * x(l,c);
                        mapa(destino_l - mod(destino_l,1) + sign(destino_l), destino_c - mod(destino_c,1) + sign(destino_c)) = mapa(destino_l - mod(destino_l,1) + sign(destino_l), destino_c - mod(destino_c,1) + sign(destino_c)) + mod(abs(destino_l),1) * mod(abs(destino_c),1);
                    end
                end
            end
        end
    end
end

for l = 1:nl,
    for c = 1:nc,
        if(mapa(l,c) > 1)
            x_(l,c) = x_(l,c)/mapa(l,c);
        else
            if(mapa(l,c) > 0)
                if(flag_espelha == 0)
                    x_(l,c) = x_(l,c);%/mapa(l,c);
                else
                    x_(l,c) = x_(l,c)*mapa(l,c) + x(l,c)*(1 - mapa(l,c));
                end
            else
                if(flag_espelha == 0)
                    x_(l,c) = 0;
                else
                    x_(l,c) = x(l,c);
                end
            end
        end
    end
end
%mapa
x = x_;

% mapa
% if(flag_espelha == 1)
%     x_ = [x_(:, nl:-1:1), x_, x_(:, nl:-1:1)];
%     x_ = [x_(nl:-1:1, :); x_; x_(nl:-1:1, :)];
%     for l = 1:nl,
%         for c = 1:nc,
%             if(mapa(l,c) < 1)
%                 temp_l = dy(l,c) - sign(dy(l,c))*mod(abs(dy(l,c)),1) + sign(dy(l,c));
%                 temp_c = dx(l,c) - sign(dx(l,c))*mod(abs(dx(l,c)),1) + sign(dx(l,c));
%                 x(l,c) = x(l,c) + (1 - mapa(l,c))*x_(l - temp_l + nl, c - temp_c + nc);
%             end
% %             if(mapa(l,c) < .9)
% %                   disp(mapa(l,c))
% %             end
%         end
%     end
% end

if(flag_lexico == 1)
    x = lexico(x);
end
