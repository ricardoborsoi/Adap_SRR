% =========================================================================
% =========================================================================
% Third example: set of natural videos, real displacements,
%                estimated motion
% 
% =========================================================================
% =========================================================================



function [frames_SR,time_reg,time_alg]=outer_super_resolve_TSR_LMS(file_LR_vid,num_frames2, nr_lr, nc_lr, upsampling_fact,fname_rec_vid)
% ===================================================
% 
% file_LR_vid - file name (will full path) of .mat file with LR video named frames_LR (an nr_lr * nc_lr * T array)
% num_frames2 - number of frames of the LR video sequence
% nr_lr - number of rows of the LR frames
% nc_lr - number of columns of the LR frames
% upsampling_fact - integer upsampling factor
% fname_rec_vid - (optional) file name of the reconstructed sequence ot be saved
% 
% frames_SR - nr_hr * nc_hr * T array with the reconstructed video sequence
% time_reg - time spent on image registration
% time_alg - time spent on image SRR
% ===================================================

% Set variables:
i=1;
parent. n_frames = num_frames2;
parent. d_factor = upsampling_fact;

parent.lr_side.nr = nr_lr;
parent.lr_side.nc = nc_lr;
parent.hr_side.nr = parent.lr_side.nr * parent.d_factor;
parent.hr_side.nc = parent.lr_side.nc * parent.d_factor;


%%%%parent.hr_side.nr = 1080;
%%%%parent.hr_side.nc = 1920;
%%%%parent.lr_side.nr = parent.hr_side.nr / parent.d_factor;
%%%%parent.lr_side.nc = parent.hr_side.nc / parent.d_factor;


% do not save video if not specified
if nargin<6
    flag2_save_video = false;
    fname_rec_vid = [];
end


% =========================================================================
% Path to folder containing the image and video files

videos_path = '~/Documents/Videos_Database_new/';
images_path = '/home/Dropbox/images';



%--------------------------------------------------------------------------
% Parameters
%--------------------------------------------------------------------------
parent. flag_compute_mse    = true;             % Flag: compute the mean square reconstruction error
parent. flag_plot_SSIM      = false;
parent. flag_compute_observable_error = false;  % Compute MSE on the LR image's subspace
parent. flag_save           = false;            % Flag: save all results (images)
parent. flag_save_video     = false;            % save the video sequences
parent. flag_boundary       = 2;                % Flag: boundary conditions for LMS - (0) zero padding (1) symmetric; (2) circular; (3) replicated
parent. flag_boundary_CCD   = 2;                % Flag: boundary conditions for decimation model - (0) zero padding (1) symmetric; (2) circular; (3) replicated
parent. flag_boundary_reg   = 2;                % Flag: boundary conditions for regularization - (0) zero padding (1) symmetric; (2) circular; (3) replicated
parent. flag_algorithm_init = 2;                % Flag: algorithm initialization - (0) zeros; (1) WGN; (2) LR image interpolation 
parent. flag_is_color       = false;            % Flag: false => grayscale, true => color videos (unimplemented)

%%%%parent. n_frames            = 100; 230; 200;              % Number of HR frames
%%%%parent. d_factor            = 2;                % Decimation factor
% parent. hr_side             = 256;              % Dimension (side) of the HR images (240)
% parent. lr_side             = parent.hr_side / parent.d_factor; % Dimension (side) of the LR images
parent. n_runs              = 1;               % Number or runs (Monte Carlo Simulations)
parent. sigma_e2            = 10;               % Additive noise variance ('e' ; y = DHx + e)

parent. im_index            = 0;                % Index of the initial image file




%%%%% different sizes for rows and columns
%%%%parent.hr_side.nr = 1080;
%%%%parent.hr_side.nc = 1920;
%%%%parent.lr_side.nr = parent.hr_side.nr / parent.d_factor;
%%%%parent.lr_side.nc = parent.hr_side.nc / parent.d_factor;




% % parent. n_runs = 3;
% % i=0;
% % i=i+1; files_vid_mat_HR{i} = '~/Documents/vids_new/DAVIS17_mat/bear_HR.mat';
% % i=i+1; files_vid_mat_HR{i} = '~/Documents/vids_new/DAVIS17_mat/bus_HR.mat';
% % i=i+1; files_vid_mat_HR{i} = '~/Documents/vids_new/DAVIS17_mat/elephant_HR.mat';
% % i=i+1; files_vid_mat_HR{i} = '~/Documents/vids_new/DAVIS17_mat/horsejump-low_HR.mat';
% % i=i+1; files_vid_mat_HR{i} = '~/Documents/vids_new/DAVIS17_mat/koala_HR.mat';
% % i=i+1; files_vid_mat_HR{i} = '~/Documents/vids_new/DAVIS17_mat/mallard-fly_HR.mat';
% % i=i+1; files_vid_mat_HR{i} = '~/Documents/vids_new/DAVIS17_mat/paragliding_HR.mat';
% % i=i+1; files_vid_mat_HR{i} = '~/Documents/vids_new/DAVIS17_mat/sheep_HR.mat';
% % i=i+1; files_vid_mat_HR{i} = '~/Documents/vids_new/DAVIS17_mat/train_HR.mat';

% % i=0;
% % i=i+1; files_vid_mat_LR{i} = '~/Documents/vids_new/DAVIS17_mat/bear_LR.mat';
% % i=i+1; files_vid_mat_LR{i} = '~/Documents/vids_new/DAVIS17_mat/bus_LR.mat';
% % i=i+1; files_vid_mat_LR{i} = '~/Documents/vids_new/DAVIS17_mat/elephant_LR.mat';
% % i=i+1; files_vid_mat_LR{i} = '~/Documents/vids_new/DAVIS17_mat/horsejump-low_LR.mat';
% % i=i+1; files_vid_mat_LR{i} = '~/Documents/vids_new/DAVIS17_mat/koala_LR.mat';
% % i=i+1; files_vid_mat_LR{i} = '~/Documents/vids_new/DAVIS17_mat/mallard-fly_LR.mat';
% % i=i+1; files_vid_mat_LR{i} = '~/Documents/vids_new/DAVIS17_mat/paragliding_LR.mat';
% % i=i+1; files_vid_mat_LR{i} = '~/Documents/vids_new/DAVIS17_mat/sheep_LR.mat';
% % i=i+1; files_vid_mat_LR{i} = '~/Documents/vids_new/DAVIS17_mat/train_LR.mat';






% ----------------------------
% Motion
% ----------------------------
parent. flag_global_motion = false;     % Flag: global motion assumption
parent. flag_motion        = 7;         % Characteristic of the global motion
                                        % 0 = random walk
                                        % 1 = random walk from file
                                        % 2 = diagonal
                                        % 3 = diagonal + random walk from file  
                                        % 4 = diagonal + random walk 
                                        % 5 = square...
                                        % 6 = stairs...
                                        % 7 = natural (not synthetic) sequence 
                                        % 8 = true SRR with natural (not synthetic) sequences wiht no HR images available
                                        
parent. erro_mov    = 0;                % WGN(0,erro_mov) added to the registration result 
parent. rw_step     = 1;                % Random walk step-size (synthetic motion between frames)
parent. rw_pdf      = 'uniform';        % Pdf of random walk
parent. rw_step2    = 1;                % Outer scaling for the translational displacement magnitude


% Include test for size of HR images and max size of natural sequences
% Fix natural sequences
try
if parent.flag_motion == 7
    parent. real_vid_seq_loader_instance = real_video_sequences_loader;
% % % % %     parent. real_vid_seq_loader_instance.vid_path = videos_path;
% % % % %     parent. n_runs = parent.real_vid_seq_loader_instance.check_files_and_params(parent);
end
catch
end





% Use one different motion vectors per MC realization
parent. flag_different_motion_per_mc_run = true;

% ---------------------------------------
% Generate a flying bird outlier
% ---------------------------------------
parent. flag_generate_flying_bird = false;
parent. t_flying_bird             = 32;
parent. delta_t_flying_bird       = 3;

% ---------------------------------------
% Registration 
% ---------------------------------------
parent. flag_registr_alg = 2;           % Flag: 0 = consider known motion
                                        %       1 = registration via DFT Registration
                                        %       2 = registration via Horn and Shunck
                                        %       3 = registration via Classic+NL
                                        %       4 = Registration via Block Matching

% -------------------------------------------------------------------------
% Other variables initialization
% -------------------------------------------------------------------------
parent. H_blur    = 1;                           % Blurring mask     
parent. H_CCD     = ones(parent.d_factor+1)/((parent.d_factor+1)^2);           % CCD mask
parent. H_reg     = fspecial('laplacian',.2);      % Regularization mask


switch parent.flag_boundary_CCD
    case 0,
        parent.flag_imfilter = 0;
    case 1,
        if(mod(max(size(parent.H_CCD)),2))
            parent.flag_imfilter = 'symmetric';
        else
            parent.flag_imfilter = 'circular';
            disp('Warning!!! Symmetric boundary conditions for decimation model are not OK implemented!!!!')
            disp('...it works only for masks with odd dimensions!!')
            disp('...circular boundary conditions was setted instead of symmetric ones!!')
        end
    case 2,
        parent.flag_imfilter = 'circular';
    case 3,
        %parent.flag_imfilter = 'replicate';
        parent.flag_imfilter = 'circular';
        disp('Warning!!! Replicated boundary conditions for decimation model are not OK implemented!!!!')
        disp('...circular boundary conditions was setted instead of symmetric ones!!')
end
if(parent.flag_boundary_reg ~= 2)
    disp('Warning!!! Regularization boundary condition not implemented!!!!')
    disp('Warning!!! Boundary condition setted up to circular!!!!')
end
parent.flag_imfilter_reg = 'circular';


% -------------------------------------------------------------------------
% Generate and store motion for all realizations
% -------------------------------------------------------------------------
parent. Motion_all_runs = zeros(parent.n_frames, 2, parent.n_runs);

if parent.flag_different_motion_per_mc_run
    for i=1:parent.n_runs
        Motion = create_displacement_vector(parent.n_frames, parent.flag_motion, parent.rw_step, parent.rw_pdf, 'Motion100_1', parent.d_factor);
        Motion = round(parent.rw_step2*(Motion));
        parent.Motion_all_runs(:, :, i) = Motion;
    end
else
    Motion = create_displacement_vector(parent.n_frames, parent.flag_motion, parent.rw_step, parent.rw_pdf, 'Motion100_1', parent.d_factor);
    Motion = round(parent.rw_step2*(Motion));
    for i=1:parent.n_runs
        parent. Motion_all_runs(:, :, i) = Motion;
    end
end

% -------------------------------------------------------------------------
% Overload Motion with saved vector for reproducibility
%%%% load Motion_saved_191115.mat





% ============================================================
% Now we initialize the algorithms with the robust parameters
% 
% ============================================================

% Initialize the evaluated algorithms
i = 0;
clear algorithms

%%%% % BICUBIC INTERP
%%%% i = i+1;
%%%% algorithms{i} = bicubic_interp_time_iter;
%%%% algorithms{i}.name = 'bic';
%%%% algorithms{i}.mode = 'omoms';%'omoms';%'spline';
%%%% algorithms{i}.color = 'k:';

%%%% % LMS
%%%% i = i+1;
%%%% algorithms{i} = lms_srr_kf_time_iter;
%%%% algorithms{i}.name = 'lms';
%%%% algorithms{i}.color = '--g';
%%%% algorithms{i}.mu_lms_kf      = 4.7;
%%%% algorithms{i}.alpha_lms_kf_T = 0;
%%%% algorithms{i}.alpha_lms_kf   = 0;
%%%% algorithms{i}.K_iter_grad    = 2;

%%%% % R-LMS
%%%% i = i+1;
%%%% algorithms{i} = lms_srr_kf_time_iter;
%%%% algorithms{i}.name = 'r_lms';
%%%% algorithms{i}.color = '--r';
%%%% algorithms{i}.mu_lms_kf      = 4.2;
%%%% algorithms{i}.alpha_lms_kf_T = 0;
%%%% algorithms{i}.alpha_lms_kf   = 40e-4;
%%%% algorithms{i}.K_iter_grad = 2;

%%%% % R-LMS-KF (old, approximated reg)
%%%% i = i+1;
%%%% algorithms{i} = lms_srr_kf_time_iter;
%%%% algorithms{i}.name = 'lms_kf_old';
%%%% algorithms{i}.color = 'k';
%%%% algorithms{i}.mu_lms_kf      = 3.4;
%%%% algorithms{i}.alpha_lms_kf_T = 0.017;
%%%% algorithms{i}.alpha_lms_kf   = 1e-4;
%%%% algorithms{i}.K_iter_grad = 2;

% NEW R-LMS-KF (exact regularization)
i = i+1;
algorithms{i} = lms_srr_new_kf_exact_reg_time_iter;
algorithms{i}.name = 'lms_kf_new';
algorithms{i}.color = 'b';
algorithms{i}.flag_convolutional_Q = true;
algorithms{i}.Q = parent.H_reg; % fspecial('laplacian',.2);
% Q = eye(), wavmtx(), etc
algorithms{i}.mu_lms_kf      = 2.2;
algorithms{i}.alpha_lms_kf_T = 16;
algorithms{i}.alpha_lms_kf   = 18e-4;
algorithms{i}.K_iter_grad = 2;



%================
% Initialize
%================

% Initialize algorithms ---------------
for i=1:length(algorithms)
    algorithms{i}.initialize(parent);
end


%%
% Run algorithm
%%%%demo_EUSIPCO_alg_loadedVids
path_and_name_LRvid = file_LR_vid;
[frames_SR,time_reg,time_alg] = inner_super_resolve(parent,algorithms,path_and_name_LRvid);


if flag2_save_video == true
    try
        save(fname_rec_vid,'frames_SR');
    catch
        disp('===============================================')
        disp('WARNING: could not save reconstructed video sequence! check file name and path')
        disp('===============================================')
    end
end

clc
disp('Time spent on image registration:')
disp(time_reg)

disp('Time spent on image SRR:')
disp(time_alg)




%%%%    mkdir(strcat(files_vid_mat_LR{mc_run},'_SR_adaptive'))
    
%%%%    for i=1:length(algorithms)
%%%%        frames_SR_adap = frames_SR_all_algs{i};
%%%%        save(strcat(files_vid_mat_LR{mc_run},'_SR_adaptive/SR_results_alg_',num2str(i),'.mat'),'frames_SR_adap')
%%%%    end    
%%%%save(strcat(files_vid_mat_LR{mc_run},'_SR_adaptive/SR_results_times.mat'),'time_algs')
    







%%%%% Plot and save results (MSE)
%%%%figure, hold on
%%%%plot(10*log10(algorithms{2}.Evv),'Color',[0 .5 0]) % LMS
%%%%plot(10*log10(algorithms{3}.Evv),'r--') % R-LMS
%%%%plot(10*log10(algorithms{5}.Evv),'b') % prop 1
%%%%plot(10*log10(algorithms{4}.Evv),'k--') % prop 2
%%%% plot(10*log10(algorithms{1}.Evv),'k:') % bic
%%%% legend('LMS','R-LMS','Proposed-1','Proposed-2','Bicubic interpolation')
%%%% xlabel('Time samples (t)')
%%%% ylabel('MSE [dB]')
%%%% 
%%%% mkdir('example4')
%%%% print('example4/mse_estim_motion','-depsc')
%%%% close all


%%%% % Plot and save results (SSIM)
%%%% figure, hold on
%%%% plot(algorithms{2}.Essim,'Color',[0 .5 0]) % LMS
%%%% plot(algorithms{3}.Essim,'r--') % R-LMS
%%%% plot(algorithms{5}.Essim,'b') % prop 1
%%%% plot(algorithms{4}.Essim,'k--') % prop 2
%%%% plot(algorithms{1}.Essim,'k:') % bic
%%%% legend('LMS','R-LMS','Proposed-1','Proposed-2','Bicubic interpolation')
%%%% xlabel('Time samples (t)')
%%%% ylabel('SSIM')
%%%% 
%%%% mkdir('example4')
%%%% print('example4/ssim_estim_motion','-depsc')
%%%% close all

%%%%% ---------------------------------------
%%%%% Convert to a struct to save
%%%%ex3_real_vids_avg.Evv_lms   = 10*log10(algorithms{2}.Evv);
%%%%ex3_real_vids_avg.Evv_r_lms = 10*log10(algorithms{3}.Evv);
%%%%ex3_real_vids_avg.Evv_prop1 = 10*log10(algorithms{5}.Evv);
%%%%ex3_real_vids_avg.Evv_prop2 = 10*log10(algorithms{4}.Evv);
%%%%ex3_real_vids_avg.Evv_bic   = 10*log10(algorithms{1}.Evv);
%%%%
%%%%ex3_real_vids_avg.Essim_lms   = algorithms{2}.Essim;
%%%%ex3_real_vids_avg.Essim_r_lms = algorithms{3}.Essim;
%%%%ex3_real_vids_avg.Essim_prop1 = algorithms{5}.Essim;
%%%%ex3_real_vids_avg.Essim_prop2 = algorithms{4}.Essim;
%%%%ex3_real_vids_avg.Essim_bic   = algorithms{1}.Essim;


%%%%% ---------------------------------------
%%%%% Save execution times
%%%%
%%%%ex3_real_vids_time.time_lms   = time_algs{2};
%%%%ex3_real_vids_time.time_r_lms = time_algs{3};
%%%%ex3_real_vids_time.time_prop1 = time_algs{5};
%%%%ex3_real_vids_time.time_prop2 = time_algs{4};
%%%%ex3_real_vids_time.time_bic   = time_algs{1};



%%%%save('example4/data_real_vids_avg')












