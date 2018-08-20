function [isbt] = isbtoeplitz(M, nl_bloco, nc_bloco)
%--------------- ajuda isbtoeplitz.m --------------------------------------
%
% function [isbt] = isbtoeplitz(M, nl_bloco, nc_bloco)
%
% Verifica se a matriz é bloco-toeplitz 
%
% Pârametros de entrada:    M           = Matriz
%                           nl_bloco    = num. de linhas do bloco
%                           nc_bloco    = num. de colunas do bloco
%
% Pârametros de saída:      isbt        = (1) é bl-toeplitz
%                                         (0) não é bl-toeplitz
%
% Guilherme Holsbach Costa 
% 07/02/2006
%--------------------------------------------------------------------------

M_toeplitz = btoeplitz(M(:, 1:nc_bloco), M(1:nl_bloco, :));

if(sum(sum(abs(M - M_toeplitz))))
    isbt = 0;
else
    isbt = 1;
end