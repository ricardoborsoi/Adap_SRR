function [Coord] = create_displacement_vector(n_mov, mov, step, pdf, file_name, d_factor)
%---- ajuda para gera_vetor_movimento.m -----------------------------------
%
%
% function [Coord] = create_displacement_vector(n_mov, mov, step, pdf, file_name, d_factor)
%
% Creates a vector containing global tranlational motion coordinates.
%
% Input vv:         n_mov   = number of displacements 
%                   mov     = characteristic of the motion
%                               0 = random walk
%                               1 = random walk from file
%                               2 = diagonal
%                               3 = diagonal + random walk from file 
%                               4 = diagonal + random walk 
%                               5 = square...
%                               6 = stairs...
%
%                   step    = (variance of the) step-size 
%                   pdf     = 'normal' or 'uniform'
%                   file_name = file containing the deterministic
%                               displacement vectors
%                   d_factor  = decimation factor considered in the SRR,
%                               used in mov 5 and 6. 
%
% Otuput vv:        Coord   = coordenadas (x,y) do movimento
%
%                   Obs.: the file containing the deterministic
%                   displacement vectors must contain a vector named Coord.
%
% Guilherme Holsbach Costa 
% 03/12/2008
%--------------------------------------------------------------------------

switch mov
    case 0,
        Coord = create_random_displacement(n_mov, step, pdf);
        %Coord = ((round(random('unif',0,1,n_mov,2))-.5)*2);
    case 1,
        load(file_name);
        Coord = Motion;
    case 2,
        Coord = step*ones(n_mov, 2);
    case 3,
        load(file_name);
        Coord = Motion;
        Coord = Coord + step*ones(n_mov, 2);
    case 4,
        Coord = ones(n_mov, 2) + create_random_displacement(n_mov, step, pdf);
    case 5,
        Coord = zeros(n_mov, 2);
        if(d_factor == 2)
            Coords = [1  0;  0  1; -1  0;  0  -1];
        else if(d_factor == 4)
                Coords = [ 1  0;  1  0;  1  0;  0  1;
                          -1  0; -1  0; -1  0;  0  1;
                           1  0;  1  0;  1  0;  0  1;
                          -1  0; -1  0; -1  0;  0 -1;
                           1  0;  1  0;  1  0;  0 -1;
                          -1  0; -1  0; -1  0;  0 -1];
            end
        end
        Coord = [0 0];
        t = 2;
        index = 1;
        while (t <= n_mov),
            if(index > max(size(Coords)))
                index = 1;
            end            
            Coord(t,:) = Coords(index,:);
            t = t + 1;
            index = index + 1;
        end
        clear Coords index i;
    case 6,
        Coord = zeros(n_mov, 2);
        t = 1;
        cont = 1;
        while (t <= n_mov),
            if(cont < d_factor)
                Coord(t,:) = [1 0];
                cont = cont + 1;
            else
                Coord(t,:) = [0 1];
                cont = 1;
            end
            t = t + 1;
        end
	case 7,
        Coord = zeros(n_mov, 2);
end

Coord(1,:) = [0 0];
