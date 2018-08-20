function [Dx, Dy] = Estimate (img1, img2, level, half_window_size)
% Name:     Lucas and Kanade motion estimation
% Input:    
%           img1, img2: Two consecutive frames, must be the same size
%           level: desired pyramid level, max is 4.
% Return:   
%           Optical flow field expressed by its x, y components.
% Note:     Computer Vision Problem Set 5
%           Programmed by Zhonghao Yang  (yangzh@cc)
%           College of Computing, Georgia Institute of Technology
%           Nov, 2001

[m, n] = size (img1);
G00 = img1; G10 = img2;
if (level>0)
    G01 = reduce (G00); G11 = reduce (G10);
end
if (level>1)
    G02 = reduce (G01); G12 = reduce (G11);
end
if (level>2)
    G03 = reduce (G02); G13 = reduce (G12);
end
if (level>3)
    G04 = reduce (G03); G14 = reduce (G13);
end
l = level;

for i=level:-1:0,
    if (l == level)
        switch (l)
        case 4, Dx = zeros (size (G04)); Dy = zeros (size (G04));
        case 3, Dx = zeros (size (G03)); Dy = zeros (size (G03));
        case 2, Dx = zeros (size (G02)); Dy = zeros (size (G02));
        case 1, Dx = zeros (size (G01)); Dy = zeros (size (G01));
        case 0, Dx = zeros (size (G00)); Dy = zeros (size (G00));
        end
    else
        Dx = expand_ (Dx); Dy = expand_ (Dy);
        Dx = Dx .* 2; Dy = Dy .* 2;
    end
    switch (l)
    case 4, 
        W = warp_ (G04, Dx, Dy); 
        [Vx, Vy] = EstimateMotion (W, G14, half_window_size);
    case 3, 
        W = warp_ (G03, Dx, Dy); 
        [Vx, Vy] = EstimateMotion (W, G13, half_window_size);
    case 2, 
        W = warp_ (G02, Dx, Dy); 
        [Vx, Vy] = EstimateMotion (W, G12, half_window_size);
    case 1, 
        W = warp_ (G01, Dx, Dy); 
        [Vx, Vy] = EstimateMotion (W, G11, half_window_size);
    case 0, 
        W = warp_ (G00, Dx, Dy); 
        [Vx, Vy] = EstimateMotion (W, G10, half_window_size);
    end
    [m, n] = size (W);
    Dx(1:m, 1:n) = Dx(1:m,1:n) + Vx; Dy(1:m, 1:n) = Dy(1:m, 1:n) + Vy;
    smooth (Dx); smooth (Dy);
    l = l - 1;
end

    
    


