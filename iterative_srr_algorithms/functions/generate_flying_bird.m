function [X] = generate_flying_bird(X, mc_run, t, Motion_all_runs, t_flying_bird, delta_t_flying_bird)
%--------------------------------------------------------------------------
% Generates an outlier consisting of N/2 sized black square on the center
% of the image
% 
% INPUT:        X - Input image
% 
% OUT:          X - Output image with black square
% 
% 
% Ricardo Borsoi
% 12/06/2015
% review 02/07/2015
%--------------------------------------------------------------------------

hr_side = length(X);

% X = min(X, [256*ones(hr_side/4,hr_side); ...
%             256*ones(hr_side/2,hr_side/4) zeros(hr_side/2,hr_side/2) 256*ones(hr_side/2,hr_side/4); ...
%             256*ones(hr_side/4,hr_side)] );

if t < t_flying_bird && t >= (t_flying_bird + delta_t_flying_bird)
    disp('WARNING: Not the right time instant (flying bird generation)!!!')
    return
end

if delta_t_flying_bird == 1
    X = min(X, [256*ones(hr_side/4,hr_side); ...
                256*ones(hr_side/2,hr_side/4) zeros(hr_side/2,hr_side/2) 256*ones(hr_side/2,hr_side/4); ...
                256*ones(hr_side/4,hr_side)] );
elseif delta_t_flying_bird > 1
%     c_offset = sum( Motion_all_runs(t_flying_bird:t-1, 2, mc_run) );
%     r_offset = sum( Motion_all_runs(t_flying_bird:t-1, 1, mc_run) );
    c_offset = sum( Motion_all_runs(t_flying_bird+1:t-1+1, 2, mc_run) );
    r_offset = sum( Motion_all_runs(t_flying_bird+1:t-1+1, 1, mc_run) );

    X = min(X, [256*ones(hr_side/4 -r_offset, hr_side); ...
                            256*ones(hr_side/2 - abs(min(hr_side/4-abs(r_offset),0)), hr_side/4 -c_offset) ...
                               zeros(hr_side/2 - abs(min(hr_side/4-abs(r_offset),0)), hr_side/2 - abs(min(hr_side/4-abs(c_offset),0)) ) ...
                            256*ones(hr_side/2 - abs(min(hr_side/4-abs(r_offset),0)), hr_side/4 +c_offset); ...
                256*ones(hr_side/4 +r_offset, hr_side)] );
    
%     X = min(X, [256*ones(hr_side/4 -r_offset, hr_side); ...
%                 256*ones(hr_side/2,  hr_side/4 +c_offset) zeros(hr_side/2,hr_side/2) 256*ones(hr_side/2, hr_side/4 -c_offset); ...
%                 256*ones(hr_side/4 +r_offset, hr_side)] );
else
    disp('Error: Delta t for flying bird must be a positive integer!!!')
end
        


