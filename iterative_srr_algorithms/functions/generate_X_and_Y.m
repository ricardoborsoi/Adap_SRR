% function [X, Y, r, c, Motion_all_runs, X_old, Y_old] = generate_X_and_Y(t, mc_run, Motion_all_runs, X, X_old, Y, Y_old, I, r, c, nr, nc, ...
%                                 flag_motion, flag_generate_flying_bird, t_flying_bird, delta_t_flying_bird, hr_side, H_blur, H_CCD, flag_imfilter, sigma_e2, d_factor)
function [X, Y, r, c, Motion_all_runs, X_old, Y_old] = generate_X_and_Y(parent, t, mc_run, I, r, c, nr, nc )
% -------------------------------------------------------------------------
% Generate HR and LR images from I and motion, for time t and mc_run
% 
% 
% 
% Ricardo Borsoi
% 18/07/2015
% updated 10/01/2017 - changed input parameters to parent struct
% -------------------------------------------------------------------------

% Set X(-1) and Y(-1) as matrices of zeros
if t == 1
    X_old = zeros( parent.hr_side );
    Y_old = zeros( parent.hr_side / parent.d_factor );
end


%----------------------------------------------------
% Generate HR images 
%----------------------------------------------------
if(parent.flag_motion < 7)
    % Create a local copy of the variable for editing
    Motion_all_runs = parent.Motion_all_runs;
    
    % Generate HR images (function of motion) ---------------------
    c = c + Motion_all_runs(t,2, mc_run);
    r = r + Motion_all_runs(t,1, mc_run);
    if((c <= 0) || (r <= 0) || (c + parent.hr_side > nc) || (r + parent.hr_side > nr))
        disp('Error... Motion out of the image limits!!!!');
        c = c - Motion_all_runs(t,2, mc_run);
        r = r - Motion_all_runs(t,1, mc_run);
        Motion_all_runs(t,1, mc_run) = 0;
        Motion_all_runs(t,2, mc_run) = 0;
    end
    if(t > 1)
        X_old = parent.X;
    end
    X = I(r:(r + parent.hr_side - 1), c:(c + parent.hr_side - 1));
    
elseif (parent.flag_motion == 7)   
    % Read HR images (natural HR image sequences) -----------------
    if(t > 1)
        X_old = parent.X;
    end
    
    Motion_all_runs = parent.Motion_all_runs;
    
%     X = parent.real_video_sequences_loader.load_next_frame(parent);
    X = double( parent.real_vid_seq_loader_instance.load_next_frame(parent) );
    
% % %     if ispc
% % % %                 X = double(imread(strcat('images\sequences\foreman\frame', num2str(t), '.png')));
% % %         X = double(imread(strcat('images\sequences\mobile\frame', num2str(t), '.png')));
% % % %                 X = double(imread(strcat('images/sequences/garden/frame', num2str(t), '.png')));
% % % %                 X = double(imread(strcat('images\sequences\bus\frame', num2str(t), '.png')));
% % % %                 X = double(imread(strcat('images\sequences\coastguard\frame', num2str(t), '.png')));
% % % %                 X = double(imread(strcat('images\sequences\container\frame', num2str(290-t), '.png')));
% % % %                 X = double(imread(strcat('images\sequences\tempete\frame', num2str(t), '.png')));
% % % %                 X = X(1:hr_side, 1:hr_side);
% % %         X = X(1:hr_side, length(X)-hr_side:length(X)-1);
% % %     else
% % %         X = double(imread(strcat('images/sequences/foreman/frame', num2str(t), '.png')));
% % % %                 X = double(imread(strcat('images/sequences/mobile/frame', num2str(t), '.png')));
% % % %                 X = double(imread(strcat('images/sequences/garden/frame', num2str(t), '.png')));
% % % %                 X = double(imread(strcat('images/sequences/bus/frame', num2str(t), '.png')));
% % % %                 X = double(imread(strcat('images/sequences/coastguard/frame', num2str(t), '.png')));
% % % %                 X = double(imread(strcat('images/sequences/container/frame', num2str(t), '.png')));
% % % %                 X = double(imread(strcat('images/sequences/tempete/frame', num2str(t), '.png')));
% % % %                 X = X(1:hr_side, 1:hr_side);
% % %         X = X(1:hr_side, length(X)-hr_side:length(X)-1);
% % %     end
end

%----------------------------------------------------
% Generate a flying bird
%----------------------------------------------------
% if flag_generate_flying_bird &&  t >= t_flying_bird && t <= (t_flying_bird + delta_t_flying_bird)
if          parent.flag_generate_flying_bird ...
        &&  t >= parent.t_flying_bird ...
        &&  t <= (parent.t_flying_bird + parent.delta_t_flying_bird - 1)
    X = generate_flying_bird(X, mc_run, t, Motion_all_runs, parent.t_flying_bird, parent.delta_t_flying_bird);
end

%----------------------------------------------------
% Generate observed (LR) images 
%----------------------------------------------------
if(t > 1)
    Y_old = parent.Y;
end

if parent.flag_motion < 8
    Temp = imfilter(X, parent.H_blur, parent.flag_imfilter, 'same');
    Temp = imfilter(Temp, parent.H_CCD, parent.flag_imfilter, 'same');
    Y = Temp(1:parent.d_factor:parent.hr_side, ...
             1:parent.d_factor:parent.hr_side);
    Y = Y + sqrt(parent.sigma_e2)*randn( parent.hr_side / parent.d_factor );
else
    disp('Warning: Real image sequences not correctly implemented!!!!')
    % Read LR images (natural LR image sequences) -----------------
    Y = (imread(strcat('images\Disk\im',num2str(21-t),'.bmp')));
%             Y = (imread(strcat('images\Disk\im',num2str(t),'.bmp')));
    Y = double(Y(1:(parent.hr_side/parent.d_factor), ...
                 1:(parent.hr_side/parent.d_factor))); % n_frames = 18, n_runs = 1
end











