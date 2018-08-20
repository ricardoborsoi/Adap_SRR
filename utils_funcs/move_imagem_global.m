function [x] = move_imagem_global(x, dx, dy, nl, nc, flag_espelha)
%--------------------------------------------------------------------------
% Executa movimento global dos pixels sobre uma imagem [x_out = G*x_in].
%
%  |1 2 3|   dx = 1   |5 6 x|
%  |4 5 6|  =======>  |8 9 x| 
%  |7 8 9|   dy = 1   |x x x|
%
% O deslocamento dx e dy representam da câmera no grid das linhas e
% colunas. Isto é: x_out = x_in - dx ; y_out = y_in - dy
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
% Obs: Funciona para deslocamentos inteiros e fracionários (13/01/05)
%
% Guilherme Holsbach Costa 
% 01/12/2004
%--------------------------------------------------------------------------

flag_lexico = 0;                    % indica se a imagem de entrada foi 
                                    % passada na forma lexicográfica
if(min(size(x)) == 1)               
    x = ilexico(x, nl, nc);         %
    flag_lexico = 1;                % indica que a imagem de entrada foi 
                                    % passada na forma lexicográfica
end

if((mod(dx,1) == 0)&(mod(dy,1) == 0))
    % ---------------------------------
    % Move para dx e dy inteiros
    % ---------------------------------
    adx = abs(dx);
    ady = abs(dy);
    temp = zeros(nl + 2*ady, nc + 2*adx);
    temp(ady+1:nl+ady, adx+1:nc+adx) = x;

    if(flag_espelha == 1)                                       % 
        temp(1:ady, :) = temp(2*ady:-1:ady+1, :);               % Espelha 
        temp(nl+ady+1:nl+2*ady, :) = temp(nl+ady:-1:nl+1, :);   % imagem em 
        temp(:, 1:adx) = temp(:, 2*adx:-1:adx+1);               % torno das 
        temp(:, nc+adx+1:nc+2*adx) = temp(:, nc+adx:-1:nc+1);   % bordas
    else
        if(flag_espelha == 2)                                   % Imagem 
            temp(1:ady,:) = temp((nl+1):(ady+nl),:);            % periódica
            temp(:,1:adx) = temp(:,(nc+1):(adx+nc));            %

            temp((ady+nl)+1:(ady+nl)+ady,:) = temp((ady+1):(2*ady),:);
            temp(:, (adx+nc)+1:(adx+nc)+adx) = temp(:, (adx+1):(2*adx));
        else
            if(flag_espelha ~= 0)
                disp('Erro no espelhamento da imagem!');
            end
        end
    end        

    x = temp(ady+dy+1:ady+dy+nl, adx+dx+1:adx+dx+nc);
else
    % ---------------------------------
    % Move para dx e/ou dy fracionários
    % ---------------------------------
    
    % Redimensona imagem (cria bordas) -------------------------
    [nl, nc] = size(x);

    adx_int = abs(dx) - mod(abs(dx),1);     % Define tamanho das bordas
    adx_frac = 0;                           % em função da parte inteira
    if(mod(dx,1) ~=1)                       % e da fracionária do passo...
        adx_frac = 1;                       %
    end                                     %
    ady_int = abs(dy) - mod(abs(dy),1);     %
    ady_frac = 0;                           %
    if(mod(dy,1) ~=1)                       %
        ady_frac = 1;                       %
    end
    
    adx = adx_int + adx_frac;
    ady = ady_int + ady_frac;
        
    temp = zeros(nl + 2*ady, nc + 2*adx);   % Redimensiona imagem
    temp(ady+1:nl+ady, adx+1:nc+adx) = x;   %
    
    % Espelha bordas -------------------------------------------

    if(flag_espelha == 1)                                       % 
        temp(1:ady, :) = temp(2*ady:-1:ady+1, :);               % Espelha 
        temp(nl+ady+1:nl+2*ady, :) = temp(nl+ady:-1:nl+1, :);   % imagem em 
        temp(:, 1:adx) = temp(:, 2*adx:-1:adx+1);               % torno das 
        temp(:, nc+adx+1:nc+2*adx) = temp(:, nc+adx:-1:nc+1);   % bordas
    else
        if(flag_espelha == 2)                                   % Imagem 
            temp(1:ady,:) = temp((nl+1):(ady+nl),:);            % periódica
            temp(:,1:adx) = temp(:,(nc+1):(adx+nc));            %

            temp((ady+nl)+1:(ady+nl)+ady,:) = temp((ady+1):(2*ady),:);
            temp(:, (adx+nc)+1:(adx+nc)+adx) = temp(:, (adx+1):(2*adx));
        else
            if(flag_espelha ~= 0)
                disp('Erro no espelhamento da imagem!');
            end
        end
    end        

    % Move parte inteira do passo -------------------------------------

    dx_int = sign(dx)*adx_int;
    dy_int = sign(dy)*ady_int;
    
    temp2 = zeros(nl+2*ady_frac, nc+2*adx_frac);
    temp2 = temp(ady+1+dy_int-ady_frac:nl+ady+dy_int+ady_frac, ...
        adx+1+dx_int-adx_frac:nc+adx+dx_int+adx_frac);

    % move parte fracionária do passo ------------------------------------
    
    mascara_x = zeros(1,2);                 % Monta máscara (2x2) de 
    mascara_y = zeros(2,1);                 % ponderação entre os pixels
    mascara_x(1) = 1 - mod(abs(dx),1);      % cobertos pelo movimento
    mascara_x(2) = mod(abs(dx),1);          % fracionário de cada pixel
    mascara_y(1) = 1 - mod(abs(dy),1);      % 
    mascara_y(2) = mod(abs(dy),1);          %
    mascara = mascara_y*mascara_x;          %
    
    delta_l = 0;                            % Define o referencial da 
    delta_c = 0;                            % máscara, de acordo com o
    if(dx >= 0)                             % sentido do mvimento (+ ou -)
        if(dy < 0)                          % e define os delta necessários
            delta_l = 1;                    %
            mascara = [0 1 ; 1 0] * mascara;%
        end                                 %
    else                                    %
        mascara = mascara * [0 1 ; 1 0];    %
        delta_c = 1;                        %
        if(dy < 0)                          %
            delta_l = 1;                    %
            mascara = [0 1 ; 1 0] * mascara;%
        end                                 %
    end 

    for l = 1:nl,                           % move parte fracionária do 
        for c = 1:nc,                       % passo
            x(l,c) = sum(sum(temp2(l+1-delta_l:l+2-delta_l,...
                c+1-delta_c:c+2-delta_c) .* mascara));
        end
    end
end

if(flag_lexico == 1)
    x = lexico(x);
end
