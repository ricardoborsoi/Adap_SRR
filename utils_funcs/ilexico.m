function I = ilexico(i, nl, nc)
%---- ajuda para ilexico.m ------------------------------------------------
%
% function I = ilexico(i, nl, nc)
%
% ILEXICO Realiza a inversa da representa�ao lexicografica de uma imagem. 
% Ex: i = lexico(I); I = ilexico(i, nl, nc);
% 
% Variaveis de entrada: i  = vetor (imagem) de entrada;
%                       nl = numero de linhas de I;
%                       nc = numero de colunas de I;
%
% Guilherme Holsbach Costa 
% 05/04/2004
% Reviewed by Ricardo Borsoi on 23/03/2015
% Reviewed by Ricardo Borsoi on 24/03/2015
%--------------------------------------------------------------------------

% 1st version:
% I = zeros(nl,nc);
% for linha = 1:nl,
%     for coluna = 1:nc,
%         I(linha,coluna) = i(coluna + (linha-1)*nc);
%     end
% end


% 2nd version:
% I = zeros(nl, nc);
% for idx = 1:nl
%     I(idx,:) = i( (nc*(idx-1)+1) : (nc*(idx-1)+nc) );
% end


% 3rd version:
I = reshape(i, nc, nl)';






