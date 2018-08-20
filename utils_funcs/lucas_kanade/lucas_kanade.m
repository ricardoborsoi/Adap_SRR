function [Dx, Dy] = lucas_kanade(img1, img2, level, half_window_size)
%--------------------------------------------------------------------------
%function [Dx, Dy] = lucas_kanade(img1, img2, level, half_window_size)
%
% Estima fluxo �ptico via m�todo de Lucas and Kanade.
%
% Par�metros de entrada:    img1             = imagem de entrada I(t-1)
%                           img2             = imagem de entrada I(t)
%                           level            = n�vel da pir�mide (max = 4)
%                           half_window_size = 
%
% Par�metros de sa�da:      [Dx Dy]  = Matrizes com vetores de velocidade 
%
% Guilherme Holsbach Costa 
% 10/10/2005
%--------------------------------------------------------------------------

%addpath lucas_kanade;

[Dx, Dy] = Estimate (img1, img2, level, half_window_size);