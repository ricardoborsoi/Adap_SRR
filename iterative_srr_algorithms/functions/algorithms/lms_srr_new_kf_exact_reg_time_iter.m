classdef lms_srr_new_kf_exact_reg_time_iter < handle
% ------------------------------------------------------------------------- 
% Performs one time iteration of the EXACT-(R)-LMS-SRR-LP algorithm
% The code uses the EXACT least perturbation regularization version 
% employing a generic regularization matrix Q
%
% The least perturbation regularization considered consists of the
% EXACT version of the proximal cost function regularization
%
% This algorithm uses matrix Q (supplied by the user), and computes
% the matrix   MTXinv = (I-Q'Q/alphaT)^-1.
% Q might either by the entire matrix (flag_convolutional_Q=false), in
% which case MTXinv is computed by direct inversion, or just a convolution
% mask (flag_convolutional_Q=true), where MTXinv is also only computed 
% and stored as a convolutional mask as well.
% 
% 
% Ricardo Borsoi
% 15/04/2016
% -------------------------------------------------------------------------

    properties
        name % e.g. 'lms', 'r_lms', 'lms_kf', etc.
        color % color for plotting
        
        % Parameters
        mu_lms_kf      = 0;
        alpha_lms_kf_T = 0;
        alpha_lms_kf   = 0;
        K_iter_grad    = 0;
        
        % If Q is the full matrix or the filter mask
        flag_convolutional_Q
        Q
        MATRIXinv
        
        X_hat     %= zeros(hr_side);
        X_hat_old %= zeros(hr_side);

%         X_hat_lms_kf_old_reg
        X_hat_old_reg
        
        Evv  %= zeros(n_frames, 1);
        Essim
        Eobs %= zeros(n_frames, 1);
        error_per_iteration %= zeros(n_frames, K_iter_lms);
        
        vidObj
        
    end
    
    
    
    
    methods
        % Initialize all porcarias particular to this algorithm
        function [] = initialize(self, parent)
            if self.flag_convolutional_Q
                % invert the matrix for the block circulant case
                % Q (and MATRIXinv) is the convolution mask (e.g. Q = H_REG)
                mascara1   = conv2(self.Q, self.Q);
                ident      = zeros(size( mascara1 ));
                ident(length(self.Q),...
                      length(self.Q)) = 1;
                  
                mascara2   = ident + (1/self.alpha_lms_kf_T) * mascara1;
                mascara2inv = ifft2( (fft2(mascara2)).^-1 );
                mascara2inv = circshift(mascara2inv, [-1 -1]);
                self.MATRIXinv = mascara2inv;
            else
                % Invert full matrix
                self.Q = full(self.Q);
                if isstruct(parent.hr_side)
                    mI = eye(parent.hr_side.nr*parent.hr_side.nc);
                else
                    mI = eye(parent.hr_side^2);
                end
                self.MATRIXinv = inv( mI + (1/self.alpha_lms_kf_T) * (self.Q'*self.Q) );
            end
        end
        
        %------------------------------------------------------------------
        % Re-initializes some variables
        function [] = reset_variables(self, parent)
            self.Evv                 = zeros(parent.n_frames, 1);
            self.Essim               = zeros(parent.n_frames, 1);
            self.Eobs                = zeros(parent.n_frames, 1);
            self.error_per_iteration = zeros(parent.n_frames, self.K_iter_grad);
            if isstruct(parent.hr_side)
                self.X_hat_old = zeros(parent.hr_side.nr, parent.hr_side.nc);
                self.X_hat     = zeros(parent.hr_side.nr, parent.hr_side.nc);
            else
                self.X_hat_old = zeros(parent.hr_side);
                self.X_hat     = zeros(parent.hr_side);
            end
        end
        
        % Reset only images -----------------
        function [] = reset_images(self, parent)
            if isstruct(parent.hr_side)
                self.X_hat_old = zeros(parent.hr_side.nr, parent.hr_side.nc);
                self.X_hat     = zeros(parent.hr_side.nr, parent.hr_side.nc);
            else
                self.X_hat_old = zeros(parent.hr_side);
                self.X_hat     = zeros(parent.hr_side);
            end
        end
        
        % Post processing, from t to t+1 -----
        function [] = post_processing(self, parent)
            % X_old => xh(t-1) receives xh(t)
            self.X_hat_old = self.X_hat;
        end
        
        
        %------------------------------------------------------------------
        function [X_hat_lms_kf] = algorithm(self, parent)
            % Register the previous estimate to use on regularization
            if parent.flag_global_motion
                self.X_hat_old_reg = move_imagem_global2(self.X_hat_old, parent.Motion_hat(parent.t_record,2), parent.Motion_hat(parent.t_record,1), parent.hr_side, parent.hr_side, parent.flag_boundary);
            else
                self.X_hat_old_reg = warp_image(self.X_hat_old, parent.Motion_hat(:,:,2), parent.Motion_hat(:,:,1), parent.flag_boundary);
            end
            
            %------------------------------------------
            % Initializes local variable for the estimation
            X_hat_lms_kf = self.X_hat;
            
            if isstruct(parent.hr_side)
                % ---------------------------------------------------------
                % ---------------------------------------------------------
                % ---------------------------------------------------------
                for k = 1:self.K_iter_grad
                    %----------------------------------------------------------------------
                    % R-LMS Lagrangian
                    %----------------------------------------------------------------------
                    Temp_lms_kf   = imfilter(X_hat_lms_kf, parent.H_blur, parent.flag_imfilter, 'same');
                    Temp_lms_kf   = imfilter(Temp_lms_kf,  parent.H_CCD,  parent.flag_imfilter, 'same');
                    Dx_hat_lms_kf = Temp_lms_kf(1:parent.d_factor:end, 1:parent.d_factor:end);

                    error_lms_kf = - (parent.Y - Dx_hat_lms_kf);

                    DtError_lms_kf = zeros( parent.hr_side.nr, parent.hr_side.nc );
                    DtError_lms_kf(1 : parent.d_factor : end, ...
                                   1 : parent.d_factor : end) = error_lms_kf;

                    DtError_lms_kf = imfilter(DtError_lms_kf(end:-1:1, end:-1:1), parent.H_CCD, parent.flag_imfilter, 'same');
                    DtError_lms_kf = imfilter(DtError_lms_kf, parent.H_blur, parent.flag_imfilter, 'same');
                    DtError_lms_kf = DtError_lms_kf(end:-1:1, end:-1:1);

                    % Thikhonov regularization---------
                    regul_lms_kf_tknv = imfilter(X_hat_lms_kf, parent.H_reg, parent.flag_imfilter_reg, 'same');
                    regul_lms_kf_tknv = imfilter(regul_lms_kf_tknv(end:-1:1, end:-1:1), parent.H_reg, parent.flag_imfilter_reg, 'same');
                    regul_lms_kf_tknv = regul_lms_kf_tknv(end:-1:1, end:-1:1);

                    % Compute the lagrangian
                    lagrangian_r_lms = 2 * DtError_lms_kf    +    2 * self.alpha_lms_kf * regul_lms_kf_tknv;


                    %----------------------------------------------------------------------
                    % 1/alphaT * Q'Q G(t) xh(t-1)
                    %----------------------------------------------------------------------

                    if self.flag_convolutional_Q == false
                        regul_lms_kf_t = (1/self.alpha_lms_kf_T) * reshape( (self.Q'*(self.Q*lexico(self.X_hat_old_reg))), parent.hr_side.nc, parent.hr_side.nr)';
                    else
                        regul_lms_kf_t = imfilter(self.X_hat_old_reg, self.Q, parent.flag_imfilter_reg, 'same');
                        regul_lms_kf_t = imfilter(regul_lms_kf_t(end:-1:1, end:-1:1), self.Q, parent.flag_imfilter_reg, 'same');
                        regul_lms_kf_t = (1/self.alpha_lms_kf_T) * regul_lms_kf_t(end:-1:1, end:-1:1);
                    end



                    %----------------------------------------------------------------------
                    % (I + 1/alphaT * Q'Q )^-1 * ( xk(t) + 1/alphaT * Q'Q G(t)xh(t-1) - mu * Nabla Lr-ms(xk(t)) )
                    %----------------------------------------------------------------------

                    Temp = X_hat_lms_kf + regul_lms_kf_t  -  self.mu_lms_kf .* lagrangian_r_lms;


                    if self.flag_convolutional_Q == false
                        Temp         = self.MATRIXinv * lexico(Temp);
                        X_hat_lms_kf = ilexico(Temp, parent.hr_side.nr, parent.hr_side.nc);
                    else
                        X_hat_lms_kf = imfilter(Temp, self.MATRIXinv, parent.flag_imfilter_reg, 'same');
                    end

                    % Compute error for k-th iteration (might remove if unecessary)
                    self.error_per_iteration(parent.t_record, k) = sum(sum(  (X_hat_lms_kf - parent.X).^2  ))  ./  (parent.hr_side.nr*parent.hr_side.nc);
                end
                
                
                
                
            else
                % ---------------------------------------------------------
                % ---------------------------------------------------------
                % ---------------------------------------------------------
                for k = 1:self.K_iter_grad
                    %----------------------------------------------------------------------
                    % R-LMS Lagrangian
                    %----------------------------------------------------------------------
                    Temp_lms_kf   = imfilter(X_hat_lms_kf, parent.H_blur, parent.flag_imfilter, 'same');
                    Temp_lms_kf   = imfilter(Temp_lms_kf,  parent.H_CCD,  parent.flag_imfilter, 'same');
                    Dx_hat_lms_kf = Temp_lms_kf(1:parent.d_factor:end, 1:parent.d_factor:end);

                    error_lms_kf = - (parent.Y - Dx_hat_lms_kf);

                    DtError_lms_kf = zeros( parent.hr_side );
                    DtError_lms_kf(1 : parent.d_factor : parent.hr_side, ...
                                   1 : parent.d_factor : parent.hr_side) = error_lms_kf;

                    DtError_lms_kf = imfilter(DtError_lms_kf(end:-1:1, end:-1:1), parent.H_CCD, parent.flag_imfilter, 'same');
                    DtError_lms_kf = imfilter(DtError_lms_kf, parent.H_blur, parent.flag_imfilter, 'same');
                    DtError_lms_kf = DtError_lms_kf(end:-1:1, end:-1:1);

                    % Thikhonov regularization---------
                    regul_lms_kf_tknv = imfilter(X_hat_lms_kf, parent.H_reg, parent.flag_imfilter_reg, 'same');
                    regul_lms_kf_tknv = imfilter(regul_lms_kf_tknv(end:-1:1, end:-1:1), parent.H_reg, parent.flag_imfilter_reg, 'same');
                    regul_lms_kf_tknv = regul_lms_kf_tknv(end:-1:1, end:-1:1);

                    % Compute the lagrangian
                    lagrangian_r_lms = 2 * DtError_lms_kf    +    2 * self.alpha_lms_kf * regul_lms_kf_tknv;


                    %----------------------------------------------------------------------
                    % 1/alphaT * Q'Q G(t) xh(t-1)
                    %----------------------------------------------------------------------

                    if self.flag_convolutional_Q == false
                        regul_lms_kf_t = (1/self.alpha_lms_kf_T) * reshape( (self.Q'*(self.Q*lexico(self.X_hat_old_reg))), parent.hr_side, parent.hr_side)';
                    else
                        regul_lms_kf_t = imfilter(self.X_hat_old_reg, self.Q, parent.flag_imfilter_reg, 'same');
                        regul_lms_kf_t = imfilter(regul_lms_kf_t(end:-1:1, end:-1:1), self.Q, parent.flag_imfilter_reg, 'same');
                        regul_lms_kf_t = (1/self.alpha_lms_kf_T) * regul_lms_kf_t(end:-1:1, end:-1:1);
                    end



                    %----------------------------------------------------------------------
                    % (I + 1/alphaT * Q'Q )^-1 * ( xk(t) + 1/alphaT * Q'Q G(t)xh(t-1) - mu * Nabla Lr-ms(xk(t)) )
                    %----------------------------------------------------------------------

                    Temp = X_hat_lms_kf + regul_lms_kf_t  -  self.mu_lms_kf .* lagrangian_r_lms;


                    if self.flag_convolutional_Q == false
                        Temp         = self.MATRIXinv * lexico(Temp);
                        X_hat_lms_kf = ilexico(Temp, parent.hr_side, parent.hr_side);
                    else
                        X_hat_lms_kf = imfilter(Temp, self.MATRIXinv, parent.flag_imfilter_reg, 'same');
                    end

                    % Compute error for k-th iteration (might remove if unecessary)
                    self.error_per_iteration(parent.t_record, k) = sum(sum(  (X_hat_lms_kf - parent.X).^2  ))  ./  (parent.hr_side^2);
                end
            end
        end
        %------------------------------------------------------------------
        
        
        %------------------------------------------------------------------
        function [] = initialize_videoObj(self)
%             self.vidObj = VideoWriter( strcat('imagens/X_hat_',self.name,'.avi'), 'Uncompressed AVI');
            self.vidObj = VideoWriter( strcat('imagens/X_hat_',self.name,'.avi'), 'Grayscale AVI');
            self.vidObj.FrameRate = 10;
            open(self.vidObj);
        end
        
        %------------------------------------------------------------------
        function [] = write_videoObj(self)
            writeVideo(self.vidObj, uint8(self.X_hat));
        end
        
        %------------------------------------------------------------------
        function [] = close_videoObj(self)
            close(self.vidObj);
        end
        
    end
    
    
end
