function [M] = bcircular(Ml1)
%--------------- ajuda bcircular.m ----------------------------------------
%
% function [M] = bcircular(Mc1)
%
% Gera uma matriz bloco-circular a partir da primeira linha de blocos.
%
% Pârametros de entrada:    Ml1     = Primeira linha de blocos
%
% Pârametros de saída:      M       = Matriz bloco circular
%
% Guilherme Holsbach Costa 
% 07/02/2006
%--------------------------------------------------------------------------

[B, N] = size(Ml1);

if(N<B)
    disp('Erro em bcircular.m! Primeira linha de blocos inconsistente...')
end

M = zeros(N,N);
M(1:B,:) = Ml1;

for l = (B+1):B:N,
    M(l:(l+B-1), :) = [M((l-B):(l-1), (N-B+1):N), M((l-B):(l-1), 1:(N-B))];
end
