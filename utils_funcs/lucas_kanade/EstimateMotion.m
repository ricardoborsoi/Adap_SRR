function [Vx, Vy] = EstimateMotion (W, G1, half_window_size)
% Name:     EstimateMotion
% Input:    
%           W: warped image (count into consideration the progress already done)
%           G1: destination image
% Return:   
%           Vx, Vy: the estimation motion in this level (x, y components)
% Note:     None

[m, n] = size (W);
Vx = zeros (size (W)); Vy = zeros (size (W));

N = zeros (2*half_window_size+1, 5);

for i = 1:m,
    l = 0;
    for j = 1-half_window_size:1+half_window_size,
        l = l + 1;
        N (l,:) = getSlice (W, G1, i, j, half_window_size);
    end
    replace = 1;    
    % use a 5x5 window
    for j = 1:n, 
        t = sum (N);
        [v, d] = eig ([t(1) t(2);t(2) t(3)]);
        namda1 = d(1,1); namda2 = d(2,2);
        if (namda1 > namda2) 
            tmp = namda1; namda1 = namda2; namda2 = tmp;
            tmp1 = v (:,1); v(:,1) = v(:,2); v(:,2) = tmp1;
        end
        if (namda2 < 0.001)
            Vx (i, j) = 0; Vy (i, j) = 0;
        elseif (namda2 > 100 * namda1)
            n2 = v(1,2) * t(4) + v(2,2) * t(5);
            Vx (i,j) = n2 * v(1,2) / namda2;
            Vy (i,j) = n2 * v(2,2) / namda2;
            % Same as evaluation in Equation 10 in paper
            %n2 = v(:,2)'*[t(4); t(5)]; d2 = v(:,2)' * v(:,2);
            %Vx (i,j) = n2 * v (1,2) / d2;
            %Vy (i,j) = n2 * v (2,2) / d2;
        else
            n1 = v(1,1) * t(4) + v(2,1) * t(5);
            n2 = v(1,2) * t(4) + v(2,2) * t(5);
            Vx (i,j) = n1 * v(1,1) / namda1 + n2 * v(1,2) / namda2;
            Vy (i,j) = n1 * v(2,1) / namda1 + n2 * v(2,2) / namda2;
            % Same as evaluation in Equation 10 in paper
            %n1 = v(:,1)' * [t(4); t(5)];n2 = v(:,2)' * [t(4); t(5)]; 
            %d1 = v(:,1)' * v(:,1); d2 = v(:,2)' * v(:,2);
            %Vx (i, j) = n1 * v(1,1) / d1 + n2 * v(1,2) / d2;
            %Vy (i, j) = n1 * v(2,1) / d1 + n2 * v(2,2) / d2;
        end
        N (replace, :) = getSlice (W, G1, i, j+half_window_size+1, half_window_size);
        replace = replace + 1;
        if (replace == 2 * half_window_size + 2) 
            replace = 1;
        end
    end
end

        
        