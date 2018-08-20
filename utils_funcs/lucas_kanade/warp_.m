function result = warp_ (img, Dx, Dy)
% Name:     warp
% Input:    img: original image
%           Dx, Dy: x, y component of optical flow field
% Output:   warped image
% Note:     None

[m, n] = size (img);

%y1 = zeros (m, n); x1 = zeros (m, n);
%y2 = zeros (m, n); x2 = zeros (m, n);

%[x1, y1] = meshgrid (1:n, 1:m);
%x1 = x1 - Dx; y1 = y1 - Dy;
%[x2, y2] = meshgrid (1:n, 1:m);
%result = griddata (x1, y1, img, x2, y2, 'linear');
[x,y] = meshgrid (1:n, 1:m);
x = x + Dx (1:m, 1:n); y = y + Dy (1:m,1:n);
for i=1:m,
    for j=1:n,
        if x(i,j)>n
            x(i,j) = n;
        end
        if x(i,j)<1
            x(i,j) = 1;
        end
        if y(i,j)>m
            y(i,j) = m;
        end
        if y(i,j)<1
            y(i,j) = 1;
        end
    end
end
result = interp2 (img, x, y, 'linear');