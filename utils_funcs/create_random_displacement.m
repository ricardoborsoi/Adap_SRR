function [Coord] = create_random_displacement(n_mov, step, pdf)
%---- help for create_random_displacement.m ------------------------------- 
%
% function [Coord] = create_random_displacement(n_mov, step, pdf)
%
% New version of "gera_vetor_movimento.m".
% Creates a vector containing an iid random global tranlational motion 
% coordinates.  
%
% Input vv  :       n_mov   = number of displacements 
%                   step    = (variance of the) step-size 
%                   pdf     = 'normal' or 'uniform'
%
%                   Obs.: 1) if uniform pdf is considered, integer step-sizes 
%                         implies in integer displacements (real step-sizes
%                         imlies in real displacements);
%                         2) normal pdfs accepts distinct variances (steps) 
%                         in x and y directions;
%                         3) the first displacement Coord(1,:) is forced to
%                         be [0, 0].
%
% Output vv :       Coord   = displacement coordinates (x,y) 
%
% Guilherme Holsbach Costa 
% 03/12/2008
% 
% Changes made on 05/08 (RABorsoi)
%--------------------------------------------------------------------------
if(size(pdf) == size('uniform'))
    if(mod(step,1) == 0)
        Coord = step*sign(rand(n_mov,2)-0.5);
%         Coord = round((2*step + 1) * rand(n_mov,2) - (2*step + 1)/2);
    else
        Coord = 2*step * rand(n_mov,2) - step;
    end
else 
    if(length(step) == 2)
        Coord = [sqrt(step(1))*randn(n_mov,1), sqrt(step(2))*randn(n_mov,1)];
    else
        Coord = sqrt(step)*randn(n_mov,2);
    end
end   
Coord(1,:) = [0 0];