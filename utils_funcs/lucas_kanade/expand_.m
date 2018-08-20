function result = expand_ (ori)   
% Name: 		expand
% Purpose: 	Expand operation in the Pyramid paper
%          	G(l,n) = EXPAND (G(l,n-1))
% Returns: 	the expansion result image
% Parameters: 
%				ori: the image to be expanded
% Notes: 	Apply 2 1D filtering instead of a huge 2D filtering to save time

[m,n] = size (ori);
mid = zeros (m, n);
m1 = m * 2; n1 = n * 2;			% m1, n1 is the size of expanding image
result = zeros (m1, n1);
w = generateFilter (0.4);		% et 1D filter using alpha = 0.4;
% along X direction
for j=1:m,
   t = zeros (1, n1);
   t(1:2:n1-1) = ori (j,1:n);   %  insert zero between consecutive pixel intensities
   tmp = conv ([ori(j,n) 0 t ori(j,1) 0], w);   % 1D filtering aftering properly handling border problem
   mid(j,1:n1) = 2 .* tmp (5:n1+4);     
   % Here instead of mirroring at the border, we use tiling technique
   % because the filter size is 5, so only 2 elements are needed on both side of the original signal
   % Also instead of downsampling, here we want to interpolate, so 0s are inserted before filtering.
   % along X direction: upsampling and Convolution
   % after this loop we get a image with double width
end
% along Y direction
for i=1:n1,
   t = zeros (1, m1);
   t(1:2:m1-1) = mid (1:m,i)';   % insert zero betwen consecutive pixel intensities along column
   tmp = conv([mid(m,i) 0 t mid(1,i) 0], w);
   result(1:m1,i) = 2 .* tmp (5:m1+4)';
   % Also only 2 elements along Y direction are needed on both border of 1D signal
   % Only take the valid part of the 1D convolution and dyadic downsampling      
   % Also instead of downsampling, here we want to interpolate, so 0s are inserted before filtering.
   % along Y direction: upsampling and convolution
   % after this we get the image with double width and height
end