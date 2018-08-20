function [Coord] = gera_vetor_movimento(n_mov, passo, pdf)
%---- ajuda para gera_vetor_movimento.m ----------------------------------- 
%
% function [Coord] = gera_vetor_movimento(n_mov, passo, pdf)
%
% Gera vetor com coordenadas de movimento global (translação). 
% Tipo de movimento: random walk com passo "passo".
%
% Parametros de entrada:    n_mov   = numero de matrizes de movimento 
%                                     (no <= n <= nf)
%                           passo   = passo do random walk
%                           pdf     = distribuição do passo 
%                                     ('normal' ou 'uniforme')
%
% Parametros de saida:      Coord   = coordenadas (x,y) do movimento
%
% Guilherme Holsbach Costa 
% 28/06/2004
%--------------------------------------------------------------------------
if(size(pdf) == size('uniforme'))
    if(mod(passo,1) == 0)
        Coord = round((2*passo + 1) * rand(n_mov,2) - (2*passo + 1)/2);
    else
        Coord = 2*passo * rand(n_mov,2) - passo;
    end
else 
    if(length(passo) > 1)
        Coord = [sqrt(passo(1))*randn(n_mov,1), sqrt(passo(2))*randn(n_mov,1)];
    else
        Coord = sqrt(passo)*randn(n_mov,2);
    end
end   
Coord(1,:) = [0 0];