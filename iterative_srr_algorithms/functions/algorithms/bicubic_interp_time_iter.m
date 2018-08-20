classdef bicubic_interp_time_iter < handle
% Performs the bicubic interpolation of the LR image sequence
% 
% 
% Ricardo Borsoi
% 15/04/2016
% -------------------------------------------------------------------------
    properties
        name % e.g. 'lms', 'r_lms', 'lms_fk', etc.
        color % color for plotting
        
        % Interpolation function ('spline' or 'omoms')
        mode = 'spline'
        
        
        X_hat     %= zeros(hr_side);
        X_hat_old %= zeros(hr_side);

%         X_hat_lms_kf_old_reg
        
        Evv  %= zeros(n_frames, 1);
        Essim
        Eobs %= zeros(n_frames, 1);

        
        vidObj %= [];
        
    end
    

    methods
        % Initialize all porcarias particular to this algorithm
        function [] = initialize(self, parent)
        end
        
        %------------------------------------------------------------------
        % Re-initializes some variables
        function [] = reset_variables(self, parent)
            self.Evv       = zeros(parent.n_frames, 1);
            self.Essim     = zeros(parent.n_frames, 1);
            self.Eobs      = zeros(parent.n_frames, 1);
            if isstruct(parent.hr_side)
                self.X_hat_old = zeros(parent.hr_side.nr, parent.hr_side.nc);
                self.X_hat     = zeros(parent.hr_side.nr, parent.hr_side.nc);
            else
                self.X_hat_old = zeros(parent.hr_side);
                self.X_hat     = zeros(parent.hr_side);
            end
        end
        
        
        % Reset only images
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
        function [X_hat_interp] = algorithm(self, parent)
            % Select interpolation function
            %----------------------------------
            if strcmp(self.mode, 'spline')
                if isstruct(parent.hr_side)
                    X_hat_interp = imresize(parent.Y, [parent.hr_side.nr parent.hr_side.nc],'bicubic');
                else
                    X_hat_interp = imresize(parent.Y, [parent.hr_side parent.hr_side],'bicubic');
                end
            
            %----------------------------------
            elseif strcmp(self.mode, 'omoms')
                % include omoms interpolator of Blu
                X_hat_interp = [];
                temp_interp = [];
                
                if isstruct(parent.hr_side)
                    % Upsample columns
                    for i=1:parent.lr_side.nc
                        temp_interp(:,i) = momssupersample(parent.Y(:,i), parent.d_factor, 3);
                    end
                    % upsample rows
                    for i=1:parent.hr_side.nr
                        X_hat_interp(i,:) = momssupersample(temp_interp(i,:)', parent.d_factor, 3)';
                    end
                    
                else
                    % Upsample columns
                    for i=1:parent.lr_side
                        temp_interp(:,i) = momssupersample(parent.Y(:,i), parent.d_factor, 3);
                    end
                    % upsample rows
                    for i=1:parent.hr_side
                        X_hat_interp(i,:) = momssupersample(temp_interp(i,:)', parent.d_factor, 3)';
                    end
                end
                
            %----------------------------------
            else
                error('Invalid bicubic interpolator specified!')
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













