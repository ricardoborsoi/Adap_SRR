function [I, r_ini, c_ini, nr, nc] = read_image_from_disk(mc_run, im_index, hr_side, n_frames, flag_motion)
% -------------------------------------------------------------------------
% Read and resize the selected image to generate the sequence
% 
% 
% 
% Ricardo Borsoi
% 21/04/2015
% -------------------------------------------------------------------------

% path???



I = (imread(strcat('im (',num2str(mc_run + im_index),').jpg')));

if(length(size(I)) == 3)
    I = rgb2gray(I);
end
[nr, nc] = size(I);

% I = imresize(I,5*[hr_side hr_side],'bicubic');
% imresize(I,(5*hr_side)/min([nr nc]),'bicubic');
% % % % % I = imresize(I,(hr_side + n_frames +200+1- 200)/min([nr nc]),'bicubic');
% I = imresize(I,(hr_side + n_frames +200+1- 200)/min([nr nc]),'nearest');
% I = imresize(I,0.49,'nearest');


if isstruct(hr_side)
    I = imresize(I,(max(hr_side.nr,hr_side.nc) + n_frames +200+1- 200)/min([nr nc]),'bicubic');
else
    I = imresize(I,(hr_side + n_frames +200+1- 200)/min([nr nc]),'bicubic');
end



switch flag_motion
    case {0,1},
    %    I = imresize(I,(5*hr_side)/min([nr nc]),'bicubic');
    case {2, 3, 6},
    %    I = imresize(I,(hr_side + n_frames + 10)/min([nr nc]),'bicubic');
    case {4}
    %    I = imresize(I,(hr_side + 2*n_frames)/min([nr nc]),'bicubic');
    case {5}
    %    I = imresize(I,(hr_side + 2*d_factor)/min([nr nc]),'bicubic');
    %    I = double(imresize(I,(hr_side + n_frames + 10)/min([nr nc]),'bicubic'));
end

I        = double(I);
[nr, nc] = size(I);


% initialize the position of the sliding window -----------------------
switch flag_motion
    case {0,1,5}
        if isstruct(hr_side)
            r_ini = round(nr/2-hr_side.nr/2);      % Initialization of the initial position
            c_ini = round(nc/2-hr_side.nc/2);      % in the image to start the random motion
        else
            r_ini = round(nr/2-hr_side/2);      % Initialization of the initial position
            c_ini = round(nc/2-hr_side/2);      % in the image to start the random motion
        end
        
    case {2,3,4,6} %{2,3,4,5,6}
        r_ini = 2;                          % Initialization of the initial position
        c_ini = 2;                          % in the image to start the random motion
end
r = r_ini;
c = c_ini;














