% supersamples an input array with an O-MOMS kernel at specified scale
% s - input sequence, must be a column vector
% scale - scale, must be an integer larger than 0
% degree - 0 to 5
%
% we assume periodicity when handling boundary values. this is to be
% changed in the next release (probably mirroring condition)

% changelog:
% 4/10/2010 Changed exp() to fft

% Copyright (c) 2010, Meng Wang
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without 
% modification, are permitted provided that the following conditions are 
% met:
% 
% Redistributions of source code must retain the above copyright notice, 
% this list of conditions and the following disclaimer.
% Redistributions in binary form must reproduce the above copyright notice,
% this list of conditions and the following disclaimer in the documentation
% and/or other materials provided with the distribution. Neither the name
% of CUHK nor the names of its contributors may be used to endorse or 
% promote products derived from this software without specific prior 
% written permission.
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS 
% IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
% THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
% PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
% CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
% EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
% PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
% PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
% LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
% NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

function [s_supersample] = momssupersample(s, scale, degree)
%% parse argument
if (nargin == 1)
	scale = 2;
	degree = 3;
	fprintf('No scale specified. Default to %d\n.', 2);
	fprintf('No degree specified. Default to %d\n.', 3);
elseif (nargin == 2)
	degree = 3;
	fprintf('No degree specified. Default to %d\n.', 3);
elseif (nargin == 3)
else
	fprintf('Invalid argument.\n');
end;

%% if n ~= 1, we interpolate the first column only
% todo: along m dimension for each column
[m, n] = size(s);
if (n ~= 1)
	s = s(:,1);
end;

%% allocate space for s_supersample
s_supersample_length = m*scale;
s_supersample = zeros(s_supersample_length, 1)';
s_indices_expand = (1:scale:s_supersample_length)';
s_indices = (1:1:m)';

s_fft = fft(s);

%% moms degree
switch (degree)
	case 0
		% nearest neighbor
		s_coef = s;
		s_coef_mirror = [s_coef; s_coef]; % we mirror the coefficients
		for k=0:1:floor(scale/2)-1
			s_supersample(s_indices_expand + k) = s_coef_mirror(s_indices);
		end;
		for k=floor(scale/2):1:scale-1
			s_supersample(s_indices_expand + k) = s_coef_mirror(s_indices + 1);
		end;
	case 1
		% linear
		s_coef = s;
		s_coef_mirror = [s_coef; s_coef];
		for k=0:1:(scale-1)
			width = k/scale;
			weight1 = omoms_at(width, degree);
			weight2 = 1 - weight1;
			s_supersample(s_indices_expand + k) = weight1*s_coef_mirror(s_indices) + weight2*s_coef_mirror(s_indices + 1);
		end;
	case 2
		phai_dtft = ...
            17/120*fft(circshift(padarray([1], [m-1], 0, 'post'), [-1])) + ...
            86/120*fft(circshift(padarray([1], [m-1], 0, 'post'), [0])) + ...
            17/120*fft(circshift(padarray([1], [m-1], 0, 'post'), [1]));
		s_coef = ifft(s_fft ./ phai_dtft);
		s_coef_mirror = [s_coef; s_coef; s_coef];
		s_coef_mirror_offset = m;
		for k=0:1:(scale - 1)
			width = k/scale;
			weight1 = omoms_at(width + 1, degree);
			weight2 = omoms_at(width, degree);
			weight3 = omoms_at(1 - width, degree);
			weight4 = 1 - weight1 - weight2 - weight3; % omoms_at(2 - width);
			s_supersample(s_indices_expand + k) = ...
				weight1*s_coef_mirror(s_indices + s_coef_mirror_offset - 1) + ...
				weight2*s_coef_mirror(s_indices + s_coef_mirror_offset) + ...
				weight3*s_coef_mirror(s_indices + s_coef_mirror_offset + 1) + ...
				weight4*s_coef_mirror(s_indices + s_coef_mirror_offset + 2);
		end;
	case 3
		phai_dtft = ...
            8/42*fft(circshift(padarray([1], [m-1], 0, 'post'), [-1])) + ...
            26/42*fft(circshift(padarray([1], [m-1], 0, 'post'), [0])) + ...
            8/42*fft(circshift(padarray([1], [m-1], 0, 'post'), [1]));
		s_coef = ifft(s_fft ./ phai_dtft);
		s_coef_mirror = [s_coef; s_coef; s_coef];
		s_coef_mirror_offset = m;
		for k=0:1:(scale - 1)
			width = k/scale;
			weight1 = omoms_at(width + 1, degree);
			weight2 = omoms_at(width, degree);
			weight3 = omoms_at(1 - width, degree);
			weight4 = 1 - weight1 - weight2 - weight3; % omoms_at(2 - width);
			s_supersample(s_indices_expand + k) = ...
				weight1*s_coef_mirror(s_indices + s_coef_mirror_offset - 1) + ...
				weight2*s_coef_mirror(s_indices + s_coef_mirror_offset) + ...
				weight3*s_coef_mirror(s_indices + s_coef_mirror_offset + 1) + ...
				weight4*s_coef_mirror(s_indices + s_coef_mirror_offset + 2);
		end;
	case 4
		phai_dtft = ...
            743/120960*fft(circshift(padarray([1], [m-1], 0, 'post'), [-2])) + ...
            25588/120960*fft(circshift(padarray([1], [m-1], 0, 'post'), [-1])) + ...
            11383/20160*fft(circshift(padarray([1], [m-1], 0, 'post'), [0])) + ...
            25588/120960*fft(circshift(padarray([1], [m-1], 0, 'post'), [1])) + ...
            743/120960*fft(circshift(padarray([1], [m-1], 0, 'post'), [2]));
		s_coef = ifft(s_fft ./ phai_dtft);
		s_coef_mirror = [s_coef; s_coef; s_coef];
		s_coef_mirror_offset = m;
		for k=0:1:(scale - 1)
			width = k/scale;
			weight1 = omoms_at(width + 2, degree);
			weight2 = omoms_at(width + 1, degree);
			weight3 = omoms_at(width, degree);
			weight4 = omoms_at(1 - width, degree);
			weight5 = omoms_at(2 - width, degree);
			weight6 = 1 - weight1 - weight2 - weight3 - weight4 - weight5; % omoms_at(3 - width, degree);
			s_supersample(s_indices_expand + k) = ...
				weight1*s_coef_mirror(s_indices + s_coef_mirror_offset - 2) + ...
				weight2*s_coef_mirror(s_indices + s_coef_mirror_offset - 1) + ...
				weight3*s_coef_mirror(s_indices + s_coef_mirror_offset) + ...
				weight4*s_coef_mirror(s_indices + s_coef_mirror_offset + 1) + ...
				weight5*s_coef_mirror(s_indices + s_coef_mirror_offset + 2) + ...
				weight6*s_coef_mirror(s_indices + s_coef_mirror_offset + 3);
		end;
	case 5
		phai_dtft = ...
            1/63360*fft(circshift(padarray([1], [m-1], 0, 'post'), [-3])) + ...
            850/63360*fft(circshift(padarray([1], [m-1], 0, 'post'), [-2])) + ...
            14351/63360*fft(circshift(padarray([1], [m-1], 0, 'post'), [-1])) + ...
            749/1440*fft(circshift(padarray([1], [m-1], 0, 'post'), [0])) + ...
            14351/63360*fft(circshift(padarray([1], [m-1], 0, 'post'), [1])) + ...
            850/63360*fft(circshift(padarray([1], [m-1], 0, 'post'), [2])) + ...
            1/63360*fft(circshift(padarray([1], [m-1], 0, 'post'), [3]));
		s_coef = ifft(s_fft ./ phai_dtft);
		s_coef_mirror = [s_coef; s_coef; s_coef];
		s_coef_mirror_offset = m;
		for k=0:1:(scale - 1)
			width = k/scale;
			weight1 = omoms_at(width + 3, degree);
			weight2 = omoms_at(width + 2, degree);
			weight3 = omoms_at(width + 1, degree);
			weight4 = omoms_at(width, degree);
			weight5 = omoms_at(1 - width, degree);
			weight6 = omoms_at(2 - width, degree);
			weight7 = omoms_at(3 - width, degree);
			weight8 = 1 - weight1 - weight2 - weight3 - weight4 - weight5 - weight6 - weight7; % omoms_at(4 - width, degree);
			s_supersample(s_indices_expand + k) = ...
				weight1*s_coef_mirror(s_indices + s_coef_mirror_offset - 3) + ...
				weight2*s_coef_mirror(s_indices + s_coef_mirror_offset - 2) + ...
				weight3*s_coef_mirror(s_indices + s_coef_mirror_offset - 1) + ...
				weight4*s_coef_mirror(s_indices + s_coef_mirror_offset) + ...
				weight5*s_coef_mirror(s_indices + s_coef_mirror_offset + 1) + ...
				weight6*s_coef_mirror(s_indices + s_coef_mirror_offset + 2) + ...
				weight7*s_coef_mirror(s_indices + s_coef_mirror_offset + 3) + ...
				weight8*s_coef_mirror(s_indices + s_coef_mirror_offset + 4);
		end;		
	otherwise
		fprintf('%d is not implemented.\n', degree);
end;