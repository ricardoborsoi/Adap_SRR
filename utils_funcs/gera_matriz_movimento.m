function [G] = gera_matriz_movimento(l, c, nl, nc)
%------------------------------------------------------------------------------------
%function [G] = gera_matriz_movimento(l, c, nl, nc)
%
% Gera matriz de movimento G(n) - usada na restauraçao.
%
% Parametros de entrada:    l, c     = coordenadas do movimento
%                           nl, nc   = dimensao dos quadros a serem gerados
%
% Parametros de saida:      G        = matriz de movimento a ser usada na restauraçao
%
%
% Guilherme Holsbach Costa 
% 27/06/2004
%------------------------------------------------------------------------------------

G = zeros(nl*nc, nl*nc);

nlm = 2*abs(l) + 1;          % Define num de linhas e colunas
ncm = 2*abs(c) + 1;          % da mascara n
Mascara = zeros(nlm, ncm);
Mascara(abs(l) + 1 + l, abs(c) + 1 + c) = 1;
G = matriz_conv(Mascara, nl, nc);
