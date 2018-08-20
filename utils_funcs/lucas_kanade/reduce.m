function result = reduce (ori)   % ori is the image to be reduced
% Name: 		reduce
% Purpose: 	Reduce operation in the Pyramid paper
% 				G(k) = REDUCE (G(k-1))
% Returns: 	the reduce image
% Parameters: 
%				ori: the image to be reduced;
% Notes:

[m,n] = size (ori);
mid = zeros (m, n);
m1 = round (m/2); n1 = round (n/2);   % m1, n1 is the size of reduced image
result = zeros (m1, n1);
w = generateFilter (0.4);				  % get 1D filter using alpha = 0.4
% along X direction
for j=1:m,
   tmp = conv([ori(j,n-1:n) ori(j,1:n) ori(j,1:2)], w);
   % Here instead of mirroring at the border, we use tiling technique
   % because the filter size is 5, so only 2 elements are needed on both 
   % side of the original signal to ensure proper operation
   mid(j,1:n1) = tmp(5:2:n+4);  
   % here mid is a temporary variable
   % Also you have to select the right pixel from the filtered result, which corresponds to image.
   % Only take the valid part of the convolution result and dyadic downsampling
   % along X direction: Convolution + downsampling
   % After this loop Mid is a image with half width
end
% along Y direction
for i=1:n1,
   tmp = conv([mid(m-1:m,i); mid(1:m,i); mid(1:2,i)]', w);
   % Also only 2 elements along Y direction are needed on both border of 1D signal
   result(1:m1,i) = tmp(5:2:m+4)';
   % Also you have to align to make sure that you've taken the right pixel
   % Only take the valid part of the 1D convolution and dyadic downsampling
   % along Y direction: Convolution + downsampling
   % after this loop result is the image with half width and half height, as expected.
end
