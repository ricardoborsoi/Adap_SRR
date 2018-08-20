function [Motion_hat] = register_images(Y, Y_old, Motion_all_runs, t, hr_side, mc_run, flag_registr_alg, flag_global_motion)
%--------------------------------------------------------------------------
% Register images Y and Y_old, returning motion field
% 
% 
% 
%--------------------------------------------------------------------------

if isstruct(hr_side)
    hr_sideR = hr_side.nr;
    hr_sideC = hr_side.nc;
else
    hr_sideR = hr_side;
    hr_sideC = hr_side;
end



switch flag_registr_alg
    % ---------------------------------------------------------------------
    case {0} % ------------------------------------------------------------
        % Known motion
        if flag_global_motion
            Motion_hat = Motion_all_runs(:, :, mc_run);
        else
            Motion_hat(:,:,1) = Motion(t,1, mc_run)*ones(hr_sideR, hr_sideC); %Dx
            Motion_hat(:,:,2) = Motion(t,2, mc_run)*ones(hr_sideR, hr_sideC); %Dy
        end
        
    % ---------------------------------------------------------------------
    case {1} % ------------------------------------------------------------
        % Fourier domain registration
        if ~flag_global_motion
            disp('Error: DSPreg only works with global motion!!!')
            error('Please change the motion flag to global_motion!')
        end
        im  = Y;
        im1 = Y_old;
        im  = imresize(im,  [hr_sideR hr_sideC], 'bicubic');
        im1 = imresize(im1, [hr_sideR hr_sideC], 'bicubic');
        [output, Greg]  = dftregistration(fft2(im1), fft2(im), 100);
        Motion_hat(t,:) = [output(3) output(4)];
        
    % ---------------------------------------------------------------------
    case {2} % ------------------------------------------------------------
        % Horn and Schunck registration
        
%         % Add function folder to path before starting
%         if t == 2
%             if ispc
%                 rmpath('funcs\flow_code_v2')
%                 addpath(genpath('funcs\hs'));
%             else
%                 rmpath('funcs/flow_code_v2')
%                 addpath(genpath('funcs/hs'));
%             end
%         end
        im  = Y;
        im1 = Y_old;
        im  = imresize(im,  [hr_sideR hr_sideC], 'bicubic');
        im1 = imresize(im1, [hr_sideR hr_sideC], 'bicubic');
%                     im  = X;
%                     im1 = X_old;
        % Horn and Schunck's #1
%                     [UV_hs] = estimate_flow_hs(im,im1); %,'lambda',0.1,'pyramid_levels',4,'pyramid_spacing',2);
%                     [UV_hs] = estimate_flow_hs(im,im1);%,'lambda',500,'pyramid_levels',2,'pyramid_spacing',1);
%         [UV_hs] = estimate_flow_hs(im,im1,'lambda',1000,'pyramid_levels',4,'pyramid_spacing',2);
        [UV_hs] = estimate_flow_interface(im,im1,'hs','lambda',1000,'pyramid_levels',4,'pyramid_spacing',2);
        Dx = UV_hs(:,:,1);
        Dy = UV_hs(:,:,2);
        if flag_global_motion
            Motion_hat(t,:) = [mean(mean(Dy)) mean(mean(Dx))];
        else
            Motion_hat(:,:,1) = Dy;
            Motion_hat(:,:,2) = Dx;
        end
    % ---------------------------------------------------------------------
    case {3} % ------------------------------------------------------------
        % NL optical flow registration
        % Add function folder to path before starting
        if t == 2
            if ispc
                rmpath('funcs\hs')
                addpath(genpath('funcs\flow_code_v2')); %'funcs\flow_code_v2\utils'
            else
                rmpath('funcs/hs')
                addpath(genpath('funcs/flow_code_v2')); %'funcs/flow_code_v2/utils'
            end
        end
        im  = Y;
        im1 = Y_old;
        im  = imresize(im,  [hr_sideR hr_sideC], 'bicubic');
        im1 = imresize(im1, [hr_sideR hr_sideC], 'bicubic');

        [UV_nl] = estimate_flow_interface(im,im1, 'classic+nl-fast');
                                                               % 'classic+nl-fast', 'classic+nl'
        Dx = UV_nl(:,:,1);
        Dy = UV_nl(:,:,2);
        if flag_global_motion
            Motion_hat(t,:) = [mean(mean(Dy)) mean(mean(Dx))];
        else
            Motion_hat(:,:,1) = Dy;
            Motion_hat(:,:,2) = Dx;
        end
    % ---------------------------------------------------------------------
    case {4} % ------------------------------------------------------------
        % Block Matching from MATLAB
        im  = Y;
        im1 = Y_old;
        im  = imresize(im,  [hr_sideR hr_sideC], 'bicubic');
        im1 = imresize(im1, [hr_sideR hr_sideC], 'bicubic');

        % Initializes objects
%         hbm = vision.BlockMatcher('ReferenceFrameSource', 'Input port', ...
%                                   'SearchMethod','Exhaustive', ...
%                                   'Overlap', [4 4], ...
%                                   'BlockSize', [5 5], ...
%                                   'MaximumDisplacement', [5 5], ...
%                                   'OutputValue', 'Horizontal and vertical components in complex form');
        
        hbm = vision.BlockMatcher('ReferenceFrameSource', 'Input port', ...
                                  'SearchMethod','Three-step', ...
                                  'Overlap', [4 4], ...
                                  'BlockSize', [5 5], ...
                                  'MaximumDisplacement', [5 5], ...
                                  'MatchCriteria', 'Mean absolute difference (MAD)', ... 'Mean square error (MSE)', ...
                                  'OutputValue', 'Horizontal and vertical components in complex form');
        
        motion = step(hbm, im, im1);
        Dx = real(motion);
        Dy = imag(motion);

        if flag_global_motion
            Motion_hat(t,:) = [mean(mean(Dy)) mean(mean(Dx))];
        else
            Motion_hat(:,:,1) = Dy;
            Motion_hat(:,:,2) = Dx;
        end
        
    % ---------------------------------------------------------------------
    case {5} % ------------------------------------------------------------
        % ...
        im  = Y;
        im1 = Y_old;
        im  = imresize(im,  [hr_sideR hr_sideC], 'bicubic');
        im1 = imresize(im1, [hr_sideR hr_sideC], 'bicubic');
        
        error('Unimplemented!')
% % %         
% % %         if flag_global_motion
% % %             Motion_hat(t,:) = [mean(mean(Dy)) mean(mean(Dx))];
% % %         else
% % %             Motion_hat(:,:,1) = Dy;
% % %             Motion_hat(:,:,2) = Dx;
% % %         end
        
    % ---------------------------------------------------------------------
    otherwise % -----------------------------------------------------------
        error('Invalid registration algorithm selected!!!')
%         if flag_global_motion
%             Motion_hat = Motion;
%         else
%             Motion_hat(:,:,1) = Dx;
%             Motion_hat(:,:,2) = Dy;
%         end
end
            
            
            
            