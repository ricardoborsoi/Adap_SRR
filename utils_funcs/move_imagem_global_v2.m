function [x] = move_imagem_global_v2(x, dx, dy, nl, nc, flag_espelha)
%--------------------------------------------------------------------------
% Executa movimento global dos pixels sobre uma imagem [x_out = G*x_in].
%
%  |1 2 3|   dx = 1   |5 6 x|
%  |4 5 6|  =======>  |8 9 x| 
%  |7 8 9|   dy = 1   |x x x|
%
% O deslocamento dx e dy representam o movimento da câmera no grid das 
% linhas e colunas. Isto é: x_out = x_in - dx ; y_out = y_in - dy
%
%
% Parâmetros de entrada:    x        = imagem de entrada (representação
%                                      lexicográfica ou não)
%                           dx, dy   = coordenadas do movimento
%                           nl, nc   = dimensão dos quadros a serem gerados
%                           flag_espelha = 0 -> considera a regiao fora da 
%                                               imagem como contendo zeros
%                                          1 -> espelha imagem ao redor das 
%                                               bordas
%                                          2 -> espelha a imagem de forma
%                                               periódica
%
% Parâmetros de saída:      x        = imagem de saída = imagem de entrada 
%                                      deslocada de dx e dy
%
% Obs: Funciona para deslocamentos inteiros e fracionários 
%
% Guilherme Holsbach Costa 
% 01/12/2009
%--------------------------------------------------------------------------

flag_lexico = 0;                    % indica se a imagem de entrada foi 
                                    % passada na forma lexicográfica
if(min(size(x)) == 1)               
    x = ilexico(x, nl, nc);         %
    flag_lexico = 1;                % indica que a imagem de entrada foi 
                                    % passada na forma lexicográfica
end

abs_dx = abs(dx);
abs_dy = abs(dy);
dx_int = sign(dx)*(abs_dx - mod(abs_dx,1)); % Calcula parte inteira e 
dy_int = sign(dy)*(abs_dy - mod(abs_dy,1)); % fracionária de dx e dy
dx_frac = sign(dx)*mod(abs_dx,1);           %
dy_frac = sign(dy)*mod(abs_dy,1);           %

[nl, nc] = size(x);

%----------------------------------
% Trata regiões de fronteira
%----------------------------------
abs_dx_int = abs_dx - mod(abs_dx,1);
abs_dy_int = abs_dy - mod(abs_dy,1);
if(mod(abs_dx,1) ~= 0)
    abs_dx = abs_dx_int + 1;
end
if(mod(abs_dy,1) ~= 0)
    abs_dy = abs_dy_int + 1;
end
temp = zeros(nl + 2*abs_dy, nc + 2*abs_dx);
temp(abs_dy+1:nl+abs_dy, abs_dx+1:nc+abs_dx) = x;

if(flag_espelha == 1)
    % Espelha imagem em torno das bordas
    temp(1:abs_dy, :) = temp(2*abs_dy:-1:abs_dy+1, :);
    temp(nl+abs_dy+1:nl+2*abs_dy, :) = temp(nl+abs_dy:-1:nl+1, :);
    temp(:, 1:abs_dx) = temp(:, 2*abs_dx:-1:abs_dx+1);
        temp(:, nc+abs_dx+1:nc+2*abs_dx) = temp(:, nc+abs_dx:-1:nc+1);
else
    if(flag_espelha == 2)
        % Imagem periódica
        temp(1:abs_dy,:) = temp((nl+1):(abs_dy+nl),:);
        temp(:,1:abs_dx) = temp(:,(nc+1):(abs_dx+nc));

        temp((abs_dy+nl)+1:(abs_dy+nl)+abs_dy,:) = temp((abs_dy+1):(2*abs_dy),:);
        temp(:, (abs_dx+nc)+1:(abs_dx+nc)+abs_dx) = temp(:, (abs_dx+1):(2*abs_dx));
    else
        if(flag_espelha ~= 0)
            disp('Erro no espelhamento da imagem!');
        end
    end
end

%----------------------------------
% Movimenta imagem
%----------------------------------
a1 = (1 - abs(dx_frac))*(1 - abs(dy_frac)); % Calcula pesos (ponderações de  
a2 = (1 - abs(dx_frac))*(abs(dy_frac));     % acordo com o deslocamento... 
a3 = (abs(dx_frac))*(1 - abs(dy_frac));     % (pixel de saída = soma 
a4 = (abs(dx_frac))*(abs(dy_frac));         % ponderada dos pixels de 
                                            % entrada
                                           
for linha = 1:nl,
    for coluna = 1:nc,
        x(linha,coluna) = a1 * temp(linha + abs_dy + dy_int                , coluna + abs_dx + dx_int                ) + ...
                          a2 * temp(linha + abs_dy + dy_int + sign(dy_frac), coluna + abs_dx + dx_int                ) + ...
                          a3 * temp(linha + abs_dy + dy_int                , coluna + abs_dx + dx_int + sign(dx_frac)) + ...
                          a4 * temp(linha + abs_dy + dy_int + sign(dy_frac), coluna + abs_dx + dx_int + sign(dx_frac));
    end
end

if(flag_lexico == 1)
    x = lexico(x);
end
