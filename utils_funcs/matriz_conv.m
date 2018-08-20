function [H] = matriz_conv(M, nl, nc, flag_espelha)
%---- ajuda para matriz_conv.m --------------------------------------------
%
% function [H] = matriz_conv(M, nl, nc, flag_espelha);
%
% Gera matriz de convolu��o H (filtragem: y = Hx).
%
% P�rametros de entrada:    M            = m�scara bidimensional
%                           nl, nc       = dimens�es da imagem a ser 
%                                          filtrada
%                           flag_espelha = 0 -> considera a regiao fora da 
%                                               imagem como contendo zeros
%                                          1 -> espelha a imagem em torno  
%                                               das bordas para tratar as 
%                                               fronteiras
%                                          2 -> espelha a imagem de forma
%                                               peri�dica
%                                          3 -> replica a �ltima linha da 
%                                               imagem
%                                          
%
% Par�metros de sa�da:      H            = matriz de convolu��o
%
%
% Guilherme Holsbach Costa 
% 14/06/2004
% Revis�o em 16/11/2008
%--------------------------------------------------------------------------
[nlm ncm] = size(M);                        % Verifica tamanho da m�scara

delta_lfm = round((nlm - 1)/2);             % Verifica tamanho das bordas
delta_cfm = round((ncm - 1)/2);             % para espelhar a imagem a ser 
if(~mod(nlm,2))                             % filtrada (v�lido para 
    delta_lim = delta_lfm - 1;              % m�scaras com lados 
else                                        % diferentes, par ou impar)...
    delta_lim = delta_lfm;                  %
end                                         %
if(~mod(ncm,2))                             %
    delta_cim = delta_cfm - 1;              %
else                                        %
    delta_cim = delta_cfm;                  %
end                                         %
                                         
% Gera modelo de imagem filtrada ------------------------------------------
X = sparse(nl + delta_lim + delta_lfm, nc + delta_cim + delta_cfm);
x = [1:(nl*nc)]';
X((1+delta_lim):(nl+delta_lim), (1+delta_cim):(nc+delta_cim)) = ...
    ilexico(x,nl,nc); %

% Espelha bordas do modelo de imagem --------------------------------------
if(flag_espelha == 1)
    X(:,1:delta_cim) = X(:,2*delta_cim:-1:delta_cim+1);
    X(1:delta_lim,:) = X(2*delta_lim:-1:delta_lim+1,:);
    X(:,delta_cim+nc+1:delta_cim+nc+delta_cfm) = X(:,delta_cim+nc:-1:delta_cim+nc-delta_cfm+1);
    X(delta_lim+nl+1:delta_lim+nl+delta_lfm,:) = X(delta_lim+nl:-1:delta_lim+nl-delta_lfm+1,:);
else
    if(flag_espelha == 2)
        X(:,1:delta_cim) = X(:,nc+1:delta_cim+nc);
        X(1:delta_lim,:) = X(nl+1:delta_lim+nl,:);
        X(:,delta_cim+nc+1:delta_cim+nc+delta_cfm) = X(:,delta_cim+1:delta_cim+delta_cfm);
        X(delta_lim+nl+1:delta_lim+nl+delta_lfm,:) = X(delta_lim+1:delta_lim+delta_lfm,:);
    else
        if(flag_espelha == 3)           
            for c = delta_cim:-1:1,
                X(:,c) = X(:, delta_cim + 1);
            end 
            for l = delta_lim:-1:1,
                X(l,:) = X(delta_lim + 1,:); 
            end  
            for c = (nc + delta_cim + 1):(nc + delta_cim + delta_cfm),
                X(:,c) = X(:, nc + delta_cim); 
            end 
            for l = (nl + delta_lim + 1):(nl + delta_lim + delta_lfm), 
                X(l,:) = X(nl + delta_lim,:); 
            end
        else
            if(flag_espelha ~= 0)
                disp('Erro no espelhamento da imagem!');
            end
        end
    end
end

% Monta matriz de convolu��o ----------------------------------------------
H = sparse(nl*nc, nl*nc);                    
linha = 0;
for l = (1 + delta_lim):(nl + delta_lim),
    for c = (1 + delta_cim):(nc + delta_cim),
       
        linha = linha + 1;
        for delta_l = -delta_lim:delta_lfm,
            for delta_c = -delta_cim:delta_cfm,
                elemento = X(l + delta_l, c + delta_c);
                if(elemento ~=0)
                    H(linha, elemento) = H(linha, elemento) + ...
                        M(delta_l + delta_lim + 1, delta_c + delta_cim + 1);
                end
            end
        end

    end
end





