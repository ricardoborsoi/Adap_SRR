function [Dx, Dy] = lucas_kanade(img1, img2, level, half_window_size)
%--------------------------------------------------------------------------
%function [Dx, Dy] = lucas_kanade(img1, img2, level, half_window_size)
%
% Estima fluxo óptico via método de Lucas and Kanade.
%
% Parâmetros de entrada:    img1             = imagem de entrada I(t-1)
%                           img2             = imagem de entrada I(t)
%                           level            = nível da pirâmide (max = 4)
%                           half_window_size = 
%
% Parâmetros de saída:      [Dx Dy]  = Matrizes com vetores de velocidade 
%
% Guilherme Holsbach Costa 
% 10/10/2005
%--------------------------------------------------------------------------

%addpath lucas_kanade;

[Dx, Dy] = Estimate (img1, img2, level, half_window_size);