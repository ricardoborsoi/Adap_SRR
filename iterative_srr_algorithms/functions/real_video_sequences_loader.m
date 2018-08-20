% Function/class to load video files

%
% IN:    vid_path - Path to the video folders


classdef real_video_sequences_loader < handle
% -------------------------------------------------------------------------
% 
% 
% 
% Note that:
%   "im_index" in the parent structure controls the offset of the first video.
%   If chosen appropriately with "mc_runs" = 1, any singe video can be
%   individually processed.
% 
% Ricardo Borsoi
% 16/12/2016
% -------------------------------------------------------------------------

    properties
        
        % Data about the video files available and to be processed
        vids
        vids_selected
        
        %
        vidreader_instance
        
        vid_path
        
        % Exclude initial/corrupted frames
        frame_offset = 1;
    end
    
    
    
    
% % % % % % %     if ispc
    
    
    methods
        
        %------------------------------------------------------------------
        function [num_runs_new] = check_files_and_params(self, parent)
            % Extracts file names from the directory and select those with
            % mkv extensions
            temp_fnames = dir(self.vid_path);
            j = 1;
            for i=1:length(temp_fnames)
                % Excludes problematic file names
                if length(temp_fnames(i).name) > 2
                    extension = temp_fnames(i).name(end-2:end);

                    if strcmp(extension,'mkv')
                        self.vids{j}.name = temp_fnames(i).name;
                        if ispc
                            self.vids{j}.full_vid_path = strcat(self.vid_path, '\', temp_fnames(i).name);
                        else
                            self.vids{j}.full_vid_path = strcat(self.vid_path, '/', temp_fnames(i).name);
                        end
                        j = j+1;
                    end
                end
            end
            
            % Check the size of the images and the number of frames
            for i=1:length(self.vids)
                v = VideoReader( self.vids{i}.full_vid_path );
                self.vids{i}.size = [v.Height v.Width];
                self.vids{i}.n_frames = v.NumberOfFrames;
            end
            
            % Check if videos comply to requirements and only select those which do
            if isstruct(parent.hr_side)
                max_HR_side = max(parent.hr_side.nr, parent.hr_side.nc);
            else
                max_HR_side = parent.hr_side;
            end
            
            j = 1;
            self.vids_selected = [];
            for i=1:length(self.vids)
                % Check for size of frames
                if min( self.vids{i}.size ) > max_HR_side
                    % Check for length of sequences
                    if self.vids{i}.n_frames > parent.n_frames
                        self.vids_selected{j} = self.vids{i};
                        j = j + 1;
                    end
                end
            end
            
            % Check if variable was created
            if sum(size(self.vids_selected)) == 0 % if exist(self.vids_selected) ~= 1
                % Return flag or something
                error('Error. No videos were loaded!')
            end
            
            % Return the # of OK videos to the parent class
            % to setup the # of MC runs evaluated
            num_runs_new = length(self.vids_selected);
            
            % If the index is different than 0, resize mc_runs
            % appropriately and truncate the video info
            if parent.im_index > 0
                num_runs_new = num_runs_new - parent.im_index;
                
% % % % % %                 % Re-attribute the video vector with reduced size
% % % % % %                 clear temp
% % % % % %                 for i=1:num_runs_new
% % % % % %                     temp{i} = self.vids_selected{parent.im_index + i};
% % % % % %                 end
% % % % % %                 self.vids_selected = [];
% % % % % %                 self.vids_selected = temp;
% % % % % % %                 self.vids_selected = self.vids_selected{parent.im_index+1:end};
            end
            
            % Truncate the number of MC runs
            num_runs_new = min(num_runs_new, parent.n_runs);
        end
        
        
        
        %------------------------------------------------------------------
        % Opens video file for the i-th MC run
        function [] = open_video_file(self, vid_idx)
            self.vidreader_instance = [];
            self.vidreader_instance = VideoReader( self.vids_selected{vid_idx}.full_vid_path );
        end
        
        
        
        %------------------------------------------------------------------
        % Load next frame for the current video file
        function [next_frame] = load_next_frame(self, parent)
            % 'readFrame' is only available after R2014b
            % next_frame = readFrame(self.vidreader_instance);
            % Workaround:
            next_frame = read(self.vidreader_instance, parent.t_record + self.frame_offset);
            
            % Convert to grayscale is necessary
            if parent.flag_is_color == false
                next_frame = rgb2gray(next_frame);
            end
            
            % Crop the image to the chosen size
%             next_frame = next_frame(1:parent.hr_side, 1:parent.hr_side);
            if isstruct(parent.hr_side)
                next_frame = next_frame(1:parent.hr_side.nr, end-parent.hr_side.nc+1:end);
            else
                next_frame = next_frame(1:parent.hr_side, end-parent.hr_side+1:end);
            end
        end
        
        
    end
end











% % % % Later print data (Nruns and Nframes)
% % % % Create method to print
% % % function verbose()
% % %     disp('============================================')
% % %     disp('=== Information for video processing ===')
% % %     disp('Number of frames:')
% % %     disp(parent.n_frames)
% % %     disp('Size of the images:')
% % %     disp(parent.hr_side)
% % %     disp('Number of MC runs:')
% % %     disp(parent.mc_runs)
% % %     disp('============================================')
% % % end



