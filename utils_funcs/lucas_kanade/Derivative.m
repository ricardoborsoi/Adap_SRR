function result = Derivative (img, x, y, direction)
% Name:     Derivative
% Input:    
%           img: input image
%           x, y: x, y component of the location
%           direction: 'x'/'y', the partial derivative direction
% Ouputt:   desired derivative
% Note:     None

[m, n] = size (img);
switch (direction)
case 'x', 
    if (x == 1)
        result = img (y, x+1) - img (y, x);
    elseif (x == n)
        result = img (y, x) - img (y, x-1);
    else
        result = 0.5 * (img (y, x+1) - img (y, x-1));
    end
case 'y', 
    if (y == 1)
        result = img (y+1, x) - img (y, x);
    elseif (y == m)
        result = img (y, x) - img (y-1, x);
    else 
        result = 0.5 * (img (y+1, x) - img (y-1, x));
    end
end
