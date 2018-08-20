function [isbc] = isbcircular(M, lado_bloco)
%--------------- ajuda isbcircular.m --------------------------------------
%
% function [isbc] = isbcircular(M, nl_bloco, nc_bloco)
%
% Verifica se a matriz é bloco-circular.
%
% Pârametros de entrada:    M           = Matriz
%                           lado_bloco  = lado do bloco
%
% Pârametros de saída:      isbc        = (1) é bl-circular
%                                         (0) não é bl-circular
%
% Guilherme Holsbach Costa 
% 07/02/2006
%--------------------------------------------------------------------------

M_circular = bcircular(M(1:lado_bloco, :));

if(sum(sum(abs(M - M_circular))))
    isbc = 0;
else
    isbc = 1;
end