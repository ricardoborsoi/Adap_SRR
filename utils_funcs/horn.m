%--------------------------------------------------------------------------
% Optical Flow Estimation
%
% Ref: B. K. P. Horn and B. G. Schunck, Determining optical flow, Artif. 
% Intell.17, 1981, 185–204.
%
%
% function [u v] = horn(I1, I2, n_iter, alpha)
%
% - u       : displacement in horizontal direction
% - v       : displacement in vertical direction
% - I1      : input image (matrix) 1
% - I2      : input image (matrix) 2
% - n_iter  : number of iterations to be performed
% - alpha   : lagrangean multiplier 
%
%
% Guilherme Holsbach Costa (holsbach@ieee.org)
% 21/02/2005
%--------------------------------------------------------------------------
function [u v] = horn(I1, I2, n_iter, alpha)
[nr nc] = size(I1);
if(size(I2) ~= size(I1))
    disp('Error! Input images with different sizes.');
else

    %----------------------------------------------------
    % variable initialization
    %----------------------------------------------------
    Ex = zeros(nr, nc);          
    Ey = zeros(nr, nc);          
    Et = zeros(nr, nc);          
    u  = zeros(nr, nc);          
    v  = zeros(nr, nc);         
    u_ = zeros(nr, nc);          
    v_ = zeros(nr, nc);         
    
    I1 = double(I1);            % Cast to double
    I2 = double(I2);            %

    %----------------------------------------------------
    % Luminance (E) partial first derivatives (in x, y and t) 
    %----------------------------------------------------
    for row = 1:nr-1,
        for col = 1:nc-1,
            Ex(row, col) = (I1(row+1, col) - I1(row, col) + ...
                            I1(row+1, col+1) - I1(row, col+1) + ...
                            I2(row+1, col) - I2(row, col) + ...
                            I2(row+1, col+1) - I2(row, col+1))/4;
    
            Ey(row, col) = (I1(row, col+1) - I1(row, col) + ...
                            I1(row+1, col+1) - I1(row+1, col) + ...
                            I2(row, col+1) - I2(row, col) + ...
                            I2(row+1, col+1) - I2(row+1, col))/4;
                       
            Et(row, col) = (I2(row, col) - I1(row, col) + ...
                            I2(row, col+1) - I1(row, col+1) + ...
                            I2(row+1, col) - I1(row+1, col) + ...
                            I2(row+1, col+1) - I1(row+1, col+1))/4;
        end
    end
    
%     Ex2 = zeros(nr,nc); Ey2 = zeros(nr,nc); Et2 = zeros(nr,nc);
%     Ex2(1:nr-1, 1:nc-1) = .25*(I1(1:nr-1, 2:nc  ) - I1(1:nr-1, 1:nc-1)+ ...
%         I1(2:nr  , 2:nc  ) - I1(2:nr  , 1:nc-1)+ ...
%         I2(1:nr-1, 2:nc  ) - I2(1:nr-1, 1:nc-1) + ...
%         I2(2:nr  , 2:nc  ) - I2(2:nr  , 1:nc-1));
%     Ey2(1:nr-1, 1:nc-1) = .25*(I1(2:nr  , 1:nc-1) - I1(1:nr-1, 1:nc-1) + ...
%         I1(2:nr  , 2:nc  ) - I1(1:nr-1, 2:nc  ) + ...
%         I2(2:nr  , 1:nc-1) - I2(1:nr-1, 1:nc-1) + ...
%         I2(2:nr  , 2:nc  ) - I2(1:nr-1, 2:nc  ));
%     Et2(1:nr-1, 1:nc-1) = .25*(I2(1:nr-1, 1:nc-1) - I1(1:nr-1, 1:nc-1) + ...
%         I2(2:nr  , 1:nc-1) - I1(2:nr  , 1:nc-1) + ...
%         I2(1:nr-1, 2:nc  ) - I1(1:nr-1, 2:nc  ) + ...
%         I2(2:nr  , 2:nc  ) - I1(2:nr  , 2:nc  ));
% 
%     disp(sum(abs(Ex(:) - Ex2(:)) + abs(Ey(:) - Ey2(:)) + abs(Et(:) - Et2(:))))
    
    %----------------------------------------------------
    % Algorithm iterations 
    %----------------------------------------------------
    for n = 1:n_iter,
        
        % Local averages computation
        for row = 2:nr-1,
            for col = 2:nc-1,
                u_(row, col) = (u(row, col-1) + u(row+1, col) + u(row, col+1) + u(row-1, col))/6 + ...
                    (u(row-1, col-1) + u(row+1, col-1) + u(row+1, col+1) + u(row-1, col+1))/12;
                v_(row, col) = (v(row, col-1) + v(row+1, col) + v(row, col+1) + v(row-1, col))/6 + ...
                    (v(row-1, col-1) + v(row+1, col-1) + v(row+1, col+1) + v(row-1, col+1))/12;
            end
        end
        
%         u_(2:nr-1, 2:nc-1) = (u(2:nr-1,1:nc-2) + u(3:nr,2:nc-1) + u(2:nr-1,3:nc) + u(1:nr-2,2:nc-1))./6 + ...
%             (u(1:nr-2,1:nc-2) + u(3:nr,1:nc-2) + u(3:nr,3:nc) + u(1:nr-2,3:nc))./12;
%         v_(2:nr-1, 2:nc-1) = (v(2:nr-1,1:nc-2) + v(3:nr,2:nc-1) + v(2:nr-1,3:nc) + v(1:nr-2,2:nc-1))./6 + ...
%             (v(1:nr-2,1:nc-2) + v(3:nr,1:nc-2) + v(3:nr,3:nc) + v(1:nr-2,3:nc))./12;

        % Minimization
        u = u_ - Ex.*(Ex.*u_ + Ey.*v_ + Et)./(alpha.^2 + Ex.^2 + Ey.^2);
        v = v_ - Ey.*(Ex.*u_ + Ey.*v_ + Et)./(alpha.^2 + Ex.^2 + Ey.^2);
    end
end 

       
        
        
        

% /****************************************************************/
% /* Compute y derivatives					*/
% /****************************************************************/
% calcIy(Ey,floatpic,t)
% float Ey[PIC_X][PIC_Y];
% float floatpic[FIVE][PIC_X][PIC_Y];
% int t;
% {
% int i,j;
% 
% printf("****** calculating Ey ******\n");
% for(i=startx;i<=endx;i+=step)
% for(j=starty;j<=endy;j+=step) 
% 	{
% 	Ey[i][j] = (floatpic[t][i][j+step] + floatpic[t][i+step][j+step] + 
% 		    floatpic[t+1][i][j+step] + floatpic[t+1][i+step][j+step])/4.0
% 		  -(floatpic[t][i][j] + floatpic[t][i+step][j] +
% 		    floatpic[t+1][i][j] + floatpic[t+1][i+step][j])/4.0;
% 	if(Ey[i][j] > BIG_MAG)
% 		{
% 		printf("Ey too large at i=%d j=%d Ey=%f\n",i,j,Ey[i][j]);
% 		exit(1);
% 		}
% 	}
% }
% 
% 
% 
% /****************************************************************/
% /* Compute average u value in neighbour about i,j		*/
% /****************************************************************/

% 
% /* Copy average values of neighbourhoods for boundaries */
% for(i=startx;i<=endx;i++)
% 	{
% 	ave[i][starty][0] = ave[i][starty+1][0];
% 	ave[i][endy][0] = ave[i][endy-1][0];
% 	ave[i][starty][1] = ave[i][starty+1][1];
% 	ave[i][endy][1] = ave[i][endy-1][1];
% 	}
% for(j=starty+1;j<=endy;j++)
% 	{
% 	ave[startx][j][0] = ave[startx+1][j][0];
% 	ave[endx][j][0] = ave[endx-1][j][0];
% 	ave[startx][j][1] = ave[startx+1][j][1];
% 	ave[endx][j][1] = ave[endx-1][j][1];
% 	}
% /* Corner Points */
% ave[startx][0][0] = ave[startx+1][1][0];
% ave[startx][0][1] = ave[startx+1][1][1];
% ave[startx][endy][0] = ave[startx+1][endy-1][0];
% ave[startx][endy][1] = ave[startx+1][endy-1][1];
% ave[endx][0][0] = ave[endx-1][1][0];
% ave[endx][0][1] = ave[endx-1][1][1];
% ave[endx][endy][0] = ave[endx-1][endy-1][0];
% ave[endx][endy][1] = ave[endx-1][endy-1][1];
% }
% 
% 
% 
% /****************************************************************/
% /* Compute u,v values						*/
% /****************************************************************/
% calc_vels(vels,vels1,Ex,Ey,Et)
% float vels[PIC_X][PIC_Y][2],vels1[PIC_X][PIC_Y][2];
% float Ex[PIC_X][PIC_Y],Ey[PIC_X][PIC_Y],Et[PIC_X][PIC_Y];
% {
% int i,j,k;
% float mag,ave[PIC_X][PIC_Y][2];
% 
% printf("****** Computing Velocity ******\n"); 
% fflush(stdout);
% vels_avg(vels1,ave);
% for(i=startx;i<=endx;i+=step)
% for(j=starty;j<=endy;j+=step) 
% 	{
% 	vels[i][j][0] = ave[i][j][0]-Ex[i][j]*
% 		(Ex[i][j]*ave[i][j][0]+Ey[i][j]*ave[i][j][1]+Et[i][j])
% 		/(alpha*alpha+Ex[i][j]*Ex[i][j]+Ey[i][j]*Ey[i][j]);
% 	vels[i][j][1] = ave[i][j][1]-Ey[i][j]*
% 		(Ex[i][j]*ave[i][j][0]+Ey[i][j]*ave[i][j][1]+Et[i][j])
% 		/(alpha*alpha+Ex[i][j]*Ex[i][j]+Ey[i][j]*Ey[i][j]);
% 	mag = sqrt(vels[i][j][0]*vels[i][j][0]+vels[i][j][1]*vels[i][j][1]);
% 	if(mag > 5.0 && FALSE) 
% 		{
% 		printf("Velocity magnitude of %f at %d %d is over 5.0\n",mag,i,j) ;
%                 }
% 	}
% }
% 
% 

