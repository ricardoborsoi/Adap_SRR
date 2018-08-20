function [image_reg] = multi_frame_warper(image, parent, idx_destination, idx_origin)
% -------------------------------------------------------------------------
% 
% "image" must be at time instant "idx_origin"
% 
% INPUT: 
% 
% OUTPUT:
% 
% Ricardo Borsoi
% 16/03/2017
% -------------------------------------------------------------------------







% When registering to the same frame, do nothing --------------------------
if     idx_destination == idx_origin
    image_reg = image;
    return;

% "forward" motion --------------------------------------------------------
elseif idx_destination > idx_origin
    % Initialize temporary variable
    temp = image;
    
    % Apply the registration recursively 
    for k=(idx_origin+1):idx_destination
        % temp = warp (G(i),temp);
        if parent.flag_global_motion
            temp = move_imagem_global2(temp, parent.Motion_hat(k,2), parent.Motion_hat(k,1), parent.hr_side, parent.hr_side, parent.flag_boundary);
        else
            temp = warp_image(temp, parent.Motion_hat{k}(:,:,2), parent.Motion_hat{k}(:,:,1), parent.flag_boundary); % warp_image(im1, Dx, Dy, flag_boundary)
        end
    end
    % Attribute the result
	image_reg = temp;

% "backward" motion -------------------------------------------------------
elseif idx_destination < idx_origin
    % Initialize temporary variable
    temp = image;
    
    % Apply the registration recursively with negative motion
    for k=idx_origin:-1:(idx_destination+1)
        % temp = warp (-G(k),temp);
        if parent.flag_global_motion
            temp = move_imagem_global2(temp, -parent.Motion_hat(k,2), -parent.Motion_hat(k,1), parent.hr_side, parent.hr_side, parent.flag_boundary);
        else
            temp = warp_image(temp, -parent.Motion_hat{k}(:,:,2), -parent.Motion_hat{k}(:,:,1), parent.flag_boundary); % warp_image(im1, Dx, Dy, flag_boundary)
        end
    end
    % Attribute the result
    image_reg = temp;
end
    
    
end



