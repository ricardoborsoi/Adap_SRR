function [X] = gera_imagens(I, Movimento, lado, l_ini, c_ini)
%---- ajuda para gera_imagens.m -------------------------------------------
%
% function [X] = gera_imagens(I, Movimento, lado, l_ini, c_ini)
%
% Gera seqüência de imagens de "alta resolução" formadas a partir da
% translação de uma janela sobre uma imagem de entrada de diemnsão maior.
%
% Parâmetros de entrada:    I            = Imagem de entrada
%                           Movimento    = coordenadas de movimento
%                           lado         = lado das imagens a serem geradas
%                           l_ini, c_ini = coordenadas iniciais da imag. de
%                                          entrada 
%
% Parâmetros de saída:      X            = conjunto de imagens geradas
%                                          (desejadas) 
%
% Guilherme Holsbach Costa 
% 28/06/2004
%--------------------------------------------------------------------------

n_imagens = length(Movimento(:,1));
X = zeros(lado, lado, n_imagens);
[nl, nc] = size(I);

l = l_ini;
c = c_ini;
for n = 1:n_imagens,
    c = c + Movimento(n,1);
    l = l + Movimento(n,2);
    if((c <= 0) | (l <= 0) | (c+lado > nc) | (l+lado > nl))
        disp('Erro... Movimento saiu da imagem!!!!');
        break;
    end
    X(:,:,n) = I(l:(l + lado - 1), c:(c + lado - 1));
end

