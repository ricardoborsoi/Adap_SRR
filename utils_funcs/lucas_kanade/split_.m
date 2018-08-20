function [odd, even] = split_ (img);
% Name:         split
% Input:        img: one frame
% Output:       
%               odd: image on the odd scanlines
%               even: image on the even scanlines
% Note:         Because the frame is taken from an interlaced video, it's advantagous to split the image
%               and handle odd and even scanlines seperately.

[m, n] = size (img);
odd = img (1:2:m, :);
even = img (2:2:m, :);
