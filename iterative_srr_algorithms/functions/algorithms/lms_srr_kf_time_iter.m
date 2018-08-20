classdef lms_srr_kf_time_iter < handle
% -------------------------------------------------------------------------
% Performs one time iteration of the (R)-LMS-SRR-(LP) algorithm
% The code may serve for the traditional LMS, thikonov regularization, and
% least perturbation versions (with and without S (S is used as default))
%
% The least perturbation regularization considered consists of the
% LINEARIZED version of the proximal cost function regularization
% 
% 
% Ricardo Borsoi
% 15/04/2016
% -------------------------------------------------------------------------
    properties
        name % e.g. 'lms', 'r_lms', 'lms_fk', etc.
        color % color for plotting
        
        % Parameters
        mu_lms_kf      = 0;
        alpha_lms_kf_T = 0;
        alpha_lms_kf   = 0;
        K_iter_grad    = 0;
        
        % Use S on temporal reg or keep equal to classical literature
        flag_use_S_reg_operator = true;
        
        
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
                % LMS KF algorithm --------------------------------------------------------
                for k = 1:self.K_iter_grad
                    Temp_lms_kf = imfilter(X_hat_lms_kf, parent.H_blur, parent.flag_imfilter, 'same');
                    Temp_lms_kf = imfilter(Temp_lms_kf,  parent.H_CCD, parent.flag_imfilter, 'same');
                    Dx_hat_lms_kf = Temp_lms_kf(1 : parent.d_factor : parent.hr_side.nr,...
                                                1 : parent.d_factor : parent.hr_side.nc);

                    error_lms_kf = parent.Y - Dx_hat_lms_kf;

                    DtError_lms_kf = zeros( parent.hr_side.nr, parent.hr_side.nc );
                    DtError_lms_kf(1 : parent.d_factor : end,...
                                   1 : parent.d_factor : end) = error_lms_kf;

                    DtError_lms_kf = imfilter(DtError_lms_kf(end:-1:1, end:-1:1), parent.H_CCD, parent.flag_imfilter, 'same');
                    DtError_lms_kf = imfilter(DtError_lms_kf, parent.H_blur, parent.flag_imfilter, 'same');
                    DtError_lms_kf = DtError_lms_kf(end:-1:1, end:-1:1);

                    DtError_lms_kf = self.mu_lms_kf.*DtError_lms_kf;
                    
                    % Regularization_lms --------------------------------------------------
                    if(self.alpha_lms_kf_T ~= 0 || self.alpha_lms_kf ~= 0)
                        regul_lms_kf      = 0;
                        
                        if self.alpha_lms_kf_T ~= 0
                            Temp_lms_regul_kf = X_hat_lms_kf - self.X_hat_old_reg;

                            %------------------------------------------------------------------
                            % Apply S'S matrix to the prior regularization term
                            if self.flag_use_S_reg_operator
                                regul_lms_kf_t = imfilter(Temp_lms_regul_kf, parent.H_reg, parent.flag_imfilter_reg, 'same');
                                regul_lms_kf_t = imfilter(regul_lms_kf_t(end:-1:1, end:-1:1), parent.H_reg, parent.flag_imfilter_reg, 'same');
                                regul_lms_kf_t = regul_lms_kf_t(end:-1:1, end:-1:1);
                                % other possible regularization options:
                                % -> [CA,CH,CV,CD] = dwt2(Temp_lms_regul_kf,'sym4');
                                %    regul_lms_kf_t = idwt2([],CH,CV,CD,'sym4', [hr_side hr_side]);
                            else
                                % Use 'S' as an identity matrix
                                regul_lms_kf_t = Temp_lms_regul_kf;
                            end %--------------------------------------------------------------

                            % Store the temporal prior value 
                            regul_lms_kf   = regul_lms_kf + self.alpha_lms_kf_T.*regul_lms_kf_t;
                        end
                        
                        % Aditional Thikhonov regularization ------------------------------
                        regul_lms_kf_tknv = imfilter(X_hat_lms_kf, parent.H_reg, parent.flag_imfilter_reg, 'same');
                        regul_lms_kf_tknv = imfilter(regul_lms_kf_tknv(end:-1:1, end:-1:1), parent.H_reg, parent.flag_imfilter_reg, 'same');
                        regul_lms_kf_tknv = regul_lms_kf_tknv(end:-1:1, end:-1:1);

                        regul_lms_kf = regul_lms_kf + self.alpha_lms_kf*regul_lms_kf_tknv;
                    else
                        regul_lms_kf = 0;
                    end

                    % Update X_hat for the k-th iteration:
                    X_hat_lms_kf = X_hat_lms_kf + DtError_lms_kf - self.mu_lms_kf*regul_lms_kf;

                    % Compute error for k-th iteration (might remove if unecessary)
                    self.error_per_iteration(parent.t_record, k) = sum(sum(  (X_hat_lms_kf - parent.X).^2  ))  ./  (parent.hr_side.nr*parent.hr_side.nc);
                end
                
                
                
            else
                % ---------------------------------------------------------
                % LMS KF algorithm --------------------------------------------------------
                for k = 1:self.K_iter_grad
                    Temp_lms_kf = imfilter(X_hat_lms_kf, parent.H_blur, parent.flag_imfilter, 'same');
                    Temp_lms_kf = imfilter(Temp_lms_kf,  parent.H_CCD, parent.flag_imfilter, 'same');
                    Dx_hat_lms_kf = Temp_lms_kf(1 : parent.d_factor : parent.hr_side,...
                                                1 : parent.d_factor : parent.hr_side);

                    error_lms_kf = parent.Y - Dx_hat_lms_kf;

                    DtError_lms_kf = zeros( parent.hr_side );
                    DtError_lms_kf(1 : parent.d_factor : end,...
                                   1 : parent.d_factor : end) = error_lms_kf;

                    DtError_lms_kf = imfilter(DtError_lms_kf(end:-1:1, end:-1:1), parent.H_CCD, parent.flag_imfilter, 'same');
                    DtError_lms_kf = imfilter(DtError_lms_kf, parent.H_blur, parent.flag_imfilter, 'same');
                    DtError_lms_kf = DtError_lms_kf(end:-1:1, end:-1:1);

                    DtError_lms_kf = self.mu_lms_kf.*DtError_lms_kf;

                    % Regularization_lms --------------------------------------------------
                    if(self.alpha_lms_kf_T ~= 0 || self.alpha_lms_kf ~= 0)
                        regul_lms_kf      = 0;
                        Temp_lms_regul_kf = X_hat_lms_kf - self.X_hat_old_reg;

                        %------------------------------------------------------------------
                        % Apply S'S matrix to the prior regularization term
                        if self.flag_use_S_reg_operator
                            regul_lms_kf_t = imfilter(Temp_lms_regul_kf, parent.H_reg, parent.flag_imfilter_reg, 'same');
                            regul_lms_kf_t = imfilter(regul_lms_kf_t(end:-1:1, end:-1:1), parent.H_reg, parent.flag_imfilter_reg, 'same');
                            regul_lms_kf_t = regul_lms_kf_t(end:-1:1, end:-1:1);
                            % other possible regularization options:
                            % -> [CA,CH,CV,CD] = dwt2(Temp_lms_regul_kf,'sym4');
                            %    regul_lms_kf_t = idwt2([],CH,CV,CD,'sym4', [hr_side hr_side]);
                        else
                            % Use 'S' as an identity matrix
                            regul_lms_kf_t = Temp_lms_regul_kf;
                        end %--------------------------------------------------------------

                        % Store the temporal prior value 
                        regul_lms_kf   = regul_lms_kf + self.alpha_lms_kf_T.*regul_lms_kf_t;

                        % Aditional Thikhonov regularization ------------------------------
                        regul_lms_kf_tknv = imfilter(X_hat_lms_kf, parent.H_reg, parent.flag_imfilter_reg, 'same');
                        regul_lms_kf_tknv = imfilter(regul_lms_kf_tknv(end:-1:1, end:-1:1), parent.H_reg, parent.flag_imfilter_reg, 'same');
                        regul_lms_kf_tknv = regul_lms_kf_tknv(end:-1:1, end:-1:1);

                        regul_lms_kf = regul_lms_kf + self.alpha_lms_kf*regul_lms_kf_tknv;
                    else
                        regul_lms_kf = 0;
                    end

                    % Update X_hat for the k-th iteration:
                    X_hat_lms_kf = X_hat_lms_kf + DtError_lms_kf - self.mu_lms_kf*regul_lms_kf;

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














