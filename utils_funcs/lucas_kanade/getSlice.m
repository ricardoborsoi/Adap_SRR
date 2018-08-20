function [N] = getSlice (W, G1, i, j, half_window_size)
% Name:     getSlice
% Input:    
%           W: image used
%           G1: destination image
%           i, j: the y and x coordinate of center of the slice concerned
% Output:   
%           N: combination of 
%               Fx * Fx: the square of spatial partial derivative along x;
%               Fx * Fy: the cross product of spatial partial derivative;
%               Fy * Fy: the square of spatial partial derivative along y;
%               Fx * Ft: the product of x spatial partial derivative and difference along t;
%               Fy * Ft: the product of y spatial partial derivative and difference along t; numerator (first term) composed of x, y components and denominator (common value)
% Note:     Here Ft is taken as W - G1, this is opposite as usual to take account of the algorithm, So now the motion estimated are actual motion, 
%           That is: from original image (x,y), x+Dx(x,y), y+Dy(x,y) will be the same object 3D point in the ending image.

N = zeros (1, 5);

[m, n] = size (W);
for y=-half_window_size:half_window_size,
    Y1 = y +i;
    if (Y1 < 1)
        Y1 = Y1 + m;
    elseif (Y1 > m)
        Y1 = Y1 - m;
    end
    X1 = j;
    if (X1 < 1)
        X1 = X1 + n;
    elseif (X1 > n)
        X1 = X1 - n;
    end
    DeriX = Derivative (G1, X1, Y1, 'x'); DeriY = Derivative (G1, X1, Y1, 'y');
    N = N + [ DeriX * DeriX, ...
        DeriX * DeriY, ...
        DeriY * DeriY, ...
        DeriX * (G1 (Y1, X1) - W (Y1, X1)), ...
        DeriY * (G1 (Y1, X1) - W (Y1, X1))];
end
