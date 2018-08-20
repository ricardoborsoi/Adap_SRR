function [x] = move_imagem_global2(x, dx, dy, nl, nc, flag_espelha)
%--------------------------------------------------------------------------
% Executa movimento global dos pixels sobre uma imagem [x_out = G*x_in].
%
%  |1 2 3|   dx = 1   |5 6 x|
%  |4 5 6|  =======>  |8 9 x| 
%  |7 8 9|   dy = 1   |x x x|
%
% O deslocamento dx e dy representam da c�mera no grid das linhas e
% colunas. Isto �: x_out = x_in - dx ; y_out = y_in - dy
%
%
% Par�metros de entrada:    x        = imagem de entrada (representa��o
%                                      lexicogr�fica ou n�o)
%                           dx, dy   = coordenadas do movimento
%                           nl, nc   = dimens�o dos quadros a serem gerados
%                           flag_espelha = 0 -> considera a regiao fora da 
%                                               imagem como contendo zeros
%                                          1 -> espelha imagem ao redor das 
%                                               bordas
%                                          2 -> espelha a imagem de forma
%                                               peri�dica
%
% Par�metros de sa�da:      x        = imagem de sa�da = imagem de entrada 
%                                      deslocada de dx e dy
%
% Obs: Funciona para deslocamentos inteiros e fracion�rios
%--------------------------------------------------------------------------

flag_lexico = 0;                    % indica se a imagem de entrada foi 
                                    % passada na forma lexicogr�fica
if(min(size(x)) == 1)               
    x = ilexico(x, nl, nc);         %
    flag_lexico = 1;                % indica que a imagem de entrada foi 
                                    % passada na forma lexicogr�fica
end


if flag_espelha == 0
    flag_imfilter_warp = 0;
elseif flag_espelha == 1
    flag_imfilter_warp = 'symmetric';
elseif flag_espelha == 2
    flag_imfilter_warp = 'circular';
else
    disp('Error: unknown boundary condition, circular BVC selected')
    flag_imfilter_warp = 'circular';
end

% Calculate kernel size and central pixel
size_max  = 2*ceil(max(abs(dy), abs(dx))) + 1;
center_px = ceil(size_max/2);

kernel = zeros(size_max);

% ---------------------------------
% Move para dx e/ou dy fracion�rios
% ---------------------------------
dx_frac = mod(abs(dx),1);
dy_frac = mod(abs(dy),1);

kernel(center_px + sign(dy)*ceil(abs(dy)),  center_px + sign(dx)*ceil(abs(dx)))  = dx_frac*dy_frac; %w2
kernel(center_px + sign(dy)*floor(abs(dy)), center_px + sign(dx)*ceil(abs(dx)))  = dx_frac*(1-dy_frac); %w1
kernel(center_px + sign(dy)*ceil(abs(dy)),  center_px + sign(dx)*floor(abs(dx))) = (1-dx_frac)*dy_frac; %w3
kernel(center_px + sign(dy)*floor(abs(dy)), center_px + sign(dx)*floor(abs(dx))) = (1-dx_frac)*(1-dy_frac); %w0

% % % % % % % % % % % % % % kernel(center_px + round(dy),  center_px + round(dx))  = 1;

% kernel


x = imfilter(x, kernel, flag_imfilter_warp, 'same');


