function [M] = btoeplitz(Mc1,Ml1)
%--------------- ajuda btoeplitz.m ----------------------------------------
%
% function [M] = btoeplitz(Mc1,Ml1)
%
% Gera uma matriz bloco-toeplitz a partir da primeira coluna de blocos e 
% da primeira linha de blocos. 
%
% P�rametros de entrada:    Mc1     = Primeira coluna de blocos
%                           Ml1     = Primeira linha de blocos
%
% P�rametros de sa�da:      M       = Matriz bloco toeplitz
%
% Guilherme Holsbach Costa 
% 07/02/2006
%--------------------------------------------------------------------------

[nl, ncb] = size(Mc1);
[nlb, nc] = size(Ml1);

if(sum(sum(abs(Mc1(1:nlb,:) - Ml1(:,1:ncb)))))
    disp('Erro em btoeplitz.m! Primeiro bloco inconsistente...')
end

M = zeros(nl,nc);
M(:,1:ncb) = Mc1;
M(1:nlb,:) = Ml1;

for l = (nlb+1):nlb:nl,
    M(l:(l+nlb-1), (ncb+1):nc) = M((l-nlb):(l-1), 1:(nc-ncb));
end






