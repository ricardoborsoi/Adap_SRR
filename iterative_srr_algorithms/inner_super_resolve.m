
function [frames_SR_all_algs,time_reg,time_algs] = inner_super_resolve(parent,algorithms,path_and_name_LRvid)
% =========================================================================
% Iterative SRR framework implementation. Only works properly for 
% one video and one algorithm
% 
% parent - structure with overall parameters
% algorithms - cell array with algorithm objects (classes)
% path_and_name_LRvid - full path to the video to be loaded, with an 
%                       nr_lr * nc_lr * T array named frames_LR
% 
%
% frames_SR_all_algs{1} - nr * nc * T array with recontructed HR frames 
% time_reg - time used by the registration algorithm
% time_algs{1} - time used by the SRR algorithm
% =========================================================================

if parent.n_runs > 1
    error('This inner function is only configured to process one video at a time!')
end
if length(algorithms) > 1
    error('This inner function is only configured to process one algorithm at a time!')
end

% path and name of file with LR video as frames_LR (nr * nc * T)
files_vid_mat_LR{1} = path_and_name_LRvid; 


% reconstructed HR images, will be initialized later
%frames_SR_all_algs{1}










% =========================================================================
% Initialize algorithms variables
% =========================================================================

for i=1:length(algorithms)
    algorithms{i}.reset_variables(parent);
end

% initialize video sequences-----------------------------------------------
if parent.flag_save_video == true
    for i=1:length(algorithms)
        algorithms{i}.initialize_videoObj;
    end
end

% Initialize simulated images----------------------------------------------
if isstruct(parent.hr_side)
    parent. X     = zeros(parent.hr_side.nr, parent.hr_side.nc);
    parent. X_old = zeros(parent.hr_side.nr, parent.hr_side.nc);
    parent. Y     = zeros(parent.lr_side.nr, parent.lr_side.nc);
    parent. Y_old = zeros(parent.lr_side.nr, parent.lr_side.nc);
else
    parent. X     = zeros(parent.hr_side);
    parent. X_old = zeros(parent.hr_side);
    parent. Y     = zeros(parent.lr_side);
    parent. Y_old = zeros(parent.lr_side);
end




% Compute algorithms and registration execution time
clear time_algs
for i=1:length(algorithms)
    time_algs{i} = 0;
end

% time_algs = 0;
time_reg  = 0;






% =========================================================================
% Monte Carlo simulations for the
% iterative SRR algorithms
% =========================================================================
for mc_run = 1:parent.n_runs
    parent.mc_run_record = mc_run;
    clc, disp(['MC run: ', num2str(parent.mc_run_record)]);
%     if ~mod(mc_run,3)
%     	input('Paused. Press enter to continue...');
%     end
    


    % &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    % load files for new videos -------------------------------------------
    %% load(files_vid_mat_HR{mc_run})
    load(files_vid_mat_LR{mc_run})
    
    nFrames_read = size(frames_LR, 3);
    
    %frames_HR = cat(3, frames_HR, 255*rand(size(frames_HR,1),size(frames_HR,2),parent.n_frames-nFrames_read));
    % frames_LR = cat(3, frames_LR, 255*rand(size(frames_LR,1),size(frames_LR,2),parent.n_frames-nFrames_read));
    if ~exist('frames_LR')
        error('Could not load "frames_LR" from file!!!')
    end
    if nFrames_read ~= parent.n_frames
        error('Inconsistent video lengths in inner and outer SR functions!')
    end
    frames_HR = zeros(parent.hr_side.nr, parent.hr_side.nc, parent.n_frames);

    clear frames_SR_all_algs
    for i=1:length(algorithms)
        frames_SR_all_algs{i} = zeros(size(frames_HR));
    end
    
    % &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    
    
    
    % Read and resize the selected image to generate the sequence ---------
    if(parent.flag_motion < 7)
        [I, r, c, nr, nc] = read_image_from_disk(mc_run, parent.im_index, parent.hr_side, parent.n_frames, parent.flag_motion);
    else
        I = []; r = 0; c = 0; nr = 0; nc = 0;
        % Load video file object inside loader class
% % % % %         parent.real_vid_seq_loader_instance.open_video_file( mc_run + parent.im_index );
    end
    
    
    % Variables initialization --------------------------------------------
    if isstruct(parent.hr_side)
        parent.X_old = zeros(parent.hr_side.nr, parent.hr_side.nc);
    else
        parent.X_old = zeros(parent.hr_side);
    end
    
    
    % Reset images
    for i=1:length(algorithms)
        algorithms{i}.reset_images(parent);
    end
    
    
    
    % Iterations on t =====================================================
    for t = 1:parent.n_frames
        parent.t_record = t;
        if(~mod(t,10))
           disp(t);
        end
        
        
        % &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        if t > nFrames_read
            break;
        end
        % &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        
        %----------------------------------------------------
        % Generate HR and LR images 
        %----------------------------------------------------
%         [parent.X, parent.Y, r, c, parent.Motion_all_runs, parent.X_old, parent.Y_old] = generate_X_and_Y(parent, t, mc_run, I, r, c, nr, nc );
        
        
        % &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        % load pre-computed video sequences
        parent.X = frames_HR(:,:,t);
        parent.Y = frames_LR(:,:,t);
        
        if t == 1
            parent.X_old = zeros(parent.hr_side.nr, parent.hr_side.nc);
            parent.Y_old = zeros(parent.lr_side.nr, parent.lr_side.nc);
        else
            parent.X_old = frames_HR(:,:,t-1);
            parent.Y_old = frames_LR(:,:,t-1);
        end
        % &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
        
        
        
        %----------------------------------------------------
        % Algorithms 
        % - search for an interpolator to initialize
        %----------------------------------------------------
        
        is_interp  = [];
        idx_interp = -1;
        for i=1:length(algorithms)
            is_interp = [is_interp   isa(algorithms{i},'bicubic_interp_time_iter')];
            if isa(algorithms{i},'bicubic_interp_time_iter')
                idx_interp = i;
            end
        end
        % If there is not an interpolator or if there are too many, 
        % use spline
        if length(is_interp) > 1  ||  length(is_interp) < 1 || idx_interp == -1
            if isstruct(parent.hr_side)
                X_hat_interp = imresize(parent.Y,[parent.hr_side.nr parent.hr_side.nc],'bicubic');
            else
                X_hat_interp = imresize(parent.Y,[parent.hr_side parent.hr_side],'bicubic');
            end
        else
            X_hat_interp = algorithms{idx_interp}.algorithm(parent);
        end
        
        
        if(t == 1)
            %----------------------------------------------------
            % Algorithm time-initialization 
            %----------------------------------------------------
            switch parent.flag_algorithm_init
                case {0} % Initialize the image with zeros
                    for i=1:length(algorithms)
                        if isstruct(parent.hr_side)
                            algorithms{i}.X_hat = zeros(parent.hr_side.nr, parent.hr_side.nc);
                        else
                            algorithms{i}.X_hat = zeros(parent.hr_side);
                        end
                    end
                    
                case {1} % Initialize the image with random noise
                    if isstruct(parent.hr_side)
                        xxx = randn(parent.hr_side.nr, parent.hr_side.nc);
                    else
                        xxx = randn(parent.hr_side);
                    end
%                     load xxx;
                    for i=1:length(algorithms)
                        algorithms{i}.X_hat = 255*(xxx-min(xxx(:)))/(max(xxx(:)) - min(xxx(:)));
                    end
                    
                case {2} % Initialize images with a bicubic interpolation
                    for i=1:length(algorithms)
                        algorithms{i}.X_hat = X_hat_interp;
                    end
            end
            
            % &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
            % store SR frames
            for i=1:length(algorithms)
                frames_SR_all_algs{i}(:,:,t) = algorithms{i}.X_hat;
            end
            % &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
            
        else
            % times -------------------------------
            tic;
            %--------------------------------------------------------------
            % Register images
            %--------------------------------------------------------------
            parent.Motion_hat = register_images(parent.Y, parent.Y_old, parent.Motion_all_runs, t, parent.hr_side, mc_run, parent.flag_registr_alg, parent.flag_global_motion);
            if parent.flag_global_motion
                parent.Motion_hat(t,:)   = parent.Motion_hat(t,:) + parent.erro_mov*randn(1,2);
            else
                parent.Motion_hat(:,:,:) = parent.Motion_hat(:,:,:) + parent.erro_mov*randn(size(parent.Motion_hat));
            end
            
            % Register last HR estimation X_hat(t-1) ----------------------
            if parent.flag_global_motion
                for i=1:length(algorithms)
                    algorithms{i}.X_hat = move_imagem_global2(algorithms{i}.X_hat, parent.Motion_hat(t,2), parent.Motion_hat(t,1), parent.hr_side, parent.hr_side, parent.flag_boundary);
                end
            else
                for i=1:length(algorithms)
                    algorithms{i}.X_hat = warp_image(algorithms{i}.X_hat, parent.Motion_hat(:,:,2), parent.Motion_hat(:,:,1), parent.flag_boundary); % warp_image(im1, Dx, Dy, flag_boundary)
                end
            end

            % times -------------------------------
            time_reg = time_reg + toc;
            
            
            
            % Registration error (noise) determination --------------------
            %[vars_reg_err] = compute_registr_err(parent.flag_global_motion, parent.flag_motion, parent.flag_boundary, Motion(t,:), Motion_hat, vars_reg_err, t, X, X_old, Y, Y_old, hr_side, lr_side);

            %==============================================================
            % Perform time iteration for all
            % algorithms in the list
            %==============================================================
            
            for i=1:length(algorithms)
                tic;
                algorithms{i}.X_hat = algorithms{i}.algorithm(parent);
                time_algs{i} = time_algs{i} + toc;
            end
            
            
            % &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
            % store SR frames
            for i=1:length(algorithms)
                frames_SR_all_algs{i}(:,:,t) = algorithms{i}.X_hat;
            end
            % &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
            
            
        %==================================================================
        end % ending "else" of "if t==1" ----------------------------------
        
        
        % Perform algorithm-specific post-processing operations -----------
        % from instant "t" to "t+1"
        for i=1:length(algorithms)
            tic;
            algorithms{i}.post_processing(parent);
            time_algs{i} = time_algs{i} + toc;
        end
	
        

        
        % Compute MSE and save results if selected ------------------------
        if (parent.flag_motion < 8) && (parent.flag_compute_mse == true)
            for i=1:length(algorithms)
                V_i = algorithms{i}.X_hat - parent.X;
                if isstruct(parent.hr_side)
                    algorithms{i}.Evv(t) = algorithms{i}.Evv(t) + sum(sum(V_i.^2)) ./ (parent.n_runs*(parent.hr_side.nr*parent.hr_side.nc));
                else
                    algorithms{i}.Evv(t) = algorithms{i}.Evv(t) + sum(sum(V_i.^2)) ./ (parent.n_runs*(parent.hr_side^2));
                end
                % Compute ssim
                algorithms{i}.Essim(t) = algorithms{i}.Essim(t) + ssim_index(parent.X, algorithms{i}.X_hat) / parent.n_runs;
                
            end
            
            if parent.flag_compute_observable_error
                for i=1:length(algorithms)
                    Temp_err_obs = imfilter(algorithms{i}.X_hat, parent.H_blur, parent.flag_imfilter, 'same');
                    Temp_err_obs = imfilter(Temp_err_obs, parent.H_CCD, parent.flag_imfilter, 'same');
                    Temp_err_obs = Temp_err_obs(1 : parent.d_factor : end, 1 : parent.d_factor : end);
                    if isstruct(parent.hr_side)
                        algorithms{i}.Eobs(t) = algorithms{i}.Eobs(t) + sum(sum( (parent.Y - Temp_err_obs).^2 )) ./ (parent.n_runs*(parent.lr_side.nr*parent.lr_side.nc));
                    else
                        algorithms{i}.Eobs(t) = algorithms{i}.Eobs(t) + sum(sum( (parent.Y - Temp_err_obs).^2 )) ./ (parent.n_runs*(parent.lr_side^2));
                    end
                end
            end
        end
        
        
        
        % Save results (images) -------------------------------------------
        if ((parent.flag_save == true)) && ( t == 100 || t == 500) && ((mc_run == 1) || (mc_run == 7) || (mc_run == 9))
            % Check if is working on a WINDOWS or UNIX machine
            if ispc % ---- WINDOWS ----
                imwrite(uint8(parent.Y), strcat('imagens\Y_t', num2str(t), 'r_', num2str(mc_run), '.png'), 'png');
                if (parent.flag_motion < 8)
                    imwrite(uint8(parent.X), strcat('imagens\X_t', num2str(t), 'r_', num2str(mc_run), '.png'), 'png');
                end
                for i=1:length(algorithms)
                    imwrite(uint8(algorithms{i}.X_hat), strcat('imagens\X_hat_', algorithms{i}.name, '_t', num2str(t), 'r_', num2str(mc_run), '.png'), 'png');
                end
                
            else % ------- UNIX ----
                imwrite(uint8(parent.Y), strcat('imagens/Y_t', num2str(t), 'r_', num2str(mc_run), '.png'), 'png');
                if (flag_motion < 8)
                    imwrite(uint8(parent.X), strcat('imagens/X_t', num2str(t), 'r_', num2str(mc_run), '.png'), 'png');
                end
                for i=1:length(algorithms)
                    imwrite(uint8(algorithms{i}.X_hat), strcat('imagens/X_hat_', algorithms{i}.name, '_t', num2str(t), 'r_', num2str(mc_run), '.png'), 'png');
                end
            end
        end
        
        
        
        
        % store video sequences--------------------------------------------
        if parent.flag_save_video == true
            for i=1:length(algorithms)
                algorithms{i}.write_videoObj;
            end
            
            % Write innovations
%             if t > 1
%                 writeVideo(vidObj_innovations, uint8(warp_image(X_old, Motion_hat(:,:,2), Motion_hat(:,:,1), flag_boundary) - X));
%             end
        end
        
    end % End time iterations ---------------------------------------------
    
    
    % &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    % save SR frames
%     for i=1:length(algorithms)
%         frames_SR_all_algs{i}(:,:,t) = algorithms{i}.X_hat;
%     end
    %%%%% mkdir(strcat(files_vid_mat_LR{mc_run},'_SR_adaptive'))
    %%%%% 
    %%%%% for i=1:length(algorithms)
    %%%%%     frames_SR_adap = frames_SR_all_algs{i};
    %%%%%     save(strcat(files_vid_mat_LR{mc_run},'_SR_adaptive/SR_results_alg_',num2str(i),'.mat'),'frames_SR_adap')
    %%%%% end
    %%%%% 
    %%%%% save(strcat(files_vid_mat_LR{mc_run},'_SR_adaptive/SR_results_times.mat'),'time_algs')
    
    
%     save(strcat(files_vid_mat_LR{mc_run},'_SR_adaptive/SR_results.mat'),'frames_SR_all_algs')
    % &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    
    
    %PSNR_lms = PSNR_lms + 10*log10(255*255*hr_side*hr_side/sum(sum((X - X_hat_lms).^2)));
end




% save video sequences ----------------------------------------------------
if parent.flag_save_video == true
    for i=1:length(algorithms)
        algorithms{i}.close_videoObj;
    end
%     % Test
%     close(vidObj_innovations);
end


% save(datestr(now))



frames_SR_all_algs = frames_SR_all_algs{1};
time_algs = time_algs{1};







