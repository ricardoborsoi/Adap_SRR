function [dx dy] = mester(I1, I2, D, sigma2_e)
%--------------------------------------------------------------------------
% function [dx dy] = mester(I1, I2, D, sigma2_e)
%
% Registro de imagens para movimento translacional e global.
%
% Ref: R. Mester, M. Hotter and Robert Bosh GmbH, Robust Displacement
% Vector Estimation Including a Statistical Error Analysis, Image
% Processing and it's Applications, Proc. 5th Intl. Conf. on, pp.168-172
% July 1995.
%
% Parâmetros de entrada:    I1       = imagem de entrada I(t)
%                           I2       = imagem de entrada I(t-1)
%                           D        = matriz de degradação 
%                                      (mod. de aquisição)
%                           sigma2_2 = pot. do ruído (mod. de aquisição)
%
% Parâmetros de saída:      [dx dy]  = vetor de deslocamento 
%
% Guilherme Holsbach Costa 
% 04/03/2005
%--------------------------------------------------------------------------
[nl nc] = size(I1);
if(sigma2_e == 0)
    sigma2_e = .1;
end
    
if(size(I2) ~= size(I1))
    disp('Erro! Imagens com tamanhos diferentes.');
else
    %addpath ..\..\LMS\funcs;
    flag_espelha = 1;
    p_e = zeros(2*D+1,2*D+1);
    
    for x = -D:D,
        for y = -D:D,
            e = lexico(I1)./255 - move_imagem_global(lexico(I2), x, y, nl, nc, flag_espelha)./255;
            p_e(x+D+1,y+D+1) = exp(-.5*e'*e/sigma2_e);
        end
    end
    Pr_total = sum(sum(p_e));
    if(Pr_total ~=0)
    Pr_d = p_e./Pr_total;
    else
    Pr_d = ones(size(p_e));
    disp('Atenção! Divisão por zero em funcs/mester.m!!!')
    end

    dx = 0 ; dy = 0;
    for x = -D:D,
        for y = -D:D,
            dx = dx + x*Pr_d(x+D+1,y+D+1);
            dy = dy + y*Pr_d(x+D+1,y+D+1);
        end
    end
end 

