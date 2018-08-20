function [H] = matriz_conv(M, nl, nc, flag_espelha);
%---- ajuda para matriz_conv.m ------------------------------------------------------
%
% function [H] = matriz_conv(M, nl, nc, flag_espelha);
%
% Gera matriz de convolucao H (filtragem: y = Hx).
%
% Parametros de entrada:    M            = mascara bidimensional
%                           nl, nc       = dimensoes da imagem a ser filtrada
%                           flag_espelha = 1 -> espelha a imagem em torno das bordas 
%                                               para tratar as fronteiras
%                                          0 -> considera a regiao fora da imagem 
%                                               como contendo zeros
%
% Parametros de saida:      H            = matriz de convolucao
%
%
% Guilherme Holsbach Costa 
% 14/06/2004
%------------------------------------------------------------------------------------
[nlm ncm] = size(M);                        % Verifica tamanho da mascara

if(~mod(nlm,2))                             % Caso a mascara tenha lado par,
    if(~mod(ncm,2))                         % coloca uma linha e/ou uma coluna a mais
        temp = zeros(nlm+1, ncm+1);         % preenchida com zeros, de forma a 
        temp(2:nlm+1, 2:ncm+1) = M;         % deixar a mascara com lado impar.
    else                                    %
        temp = zeros(nlm+1, ncm);
        temp(2:nlm+1, 1:ncm) = M;
    end
    M = temp;
    [nlm ncm] = size(M);
else
    if(~mod(ncm,2))
        temp = zeros(nlm, ncm+1);
        temp(1:nlm, 2:ncm+1) = M;
        M = temp;
        [nlm ncm] = size(M);
    end
end

delta_lm = (nlm - 1)/2;                     % Verifica tamanho das bordas para espelhar a 
delta_cm = (ncm - 1)/2;                     % imagem a ser filtrada

X = zeros(nl + 2*delta_lm, nc + 2*delta_cm);% Gera modelo de imagem filtrada
x = [1:(nl*nc)]';                           %

X((1+delta_lm):(nl+delta_lm), (1+delta_cm):(nc+delta_cm)) = ilexico(x,nl,nc); %

if(flag_espelha == 1)                                       % Espelha bordas do modelo de imagem
    for c = delta_cm:-1:1,                                  %
        X(:,c) = X(:, delta_cm + 1 + (delta_cm + 1) - c);   %
    end                                                     %
    for l = delta_lm:-1:1,                                  %
        X(l,:) = X(delta_lm + 1 + (delta_lm + 1) - l,:);    %
    end                                                     %
    for c = (nc + delta_cm + 1):(nc + delta_cm + delta_cm), %
        X(:,c) = X(:, nc + delta_cm - (c - nc - delta_cm)); %
    end                                                     %
    for l = (nl + delta_lm + 1):(nl + delta_lm + delta_lm), %
        X(l,:) = X(nl + delta_lm - (l - nl - delta_lm),:);  %
    end                                                     %
end

H = sparse(nl*nc, nl*nc);                    % Monta matriz de convolucao
linha = 0;
for l = (1 + delta_lm):(nl + delta_lm),
    for c = (1 + delta_cm):(nc + delta_cm),
       
        linha = linha + 1;
        for delta_l = -delta_lm:delta_lm,
            for delta_c = -delta_cm:delta_cm,
                elemento = X(l + delta_l, c + delta_c);
                if(elemento ~=0)
                H(linha, elemento) = H(linha, elemento) + ...
                    M(delta_l + delta_lm + 1, delta_c + delta_cm + 1);
                end
            end
        end

    end
end
