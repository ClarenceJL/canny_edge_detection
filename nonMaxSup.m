function edgeMap = nonMaxSup( J, theta )
%% nonMaxSup
%   Non maximum supression
%  Input:
%   J - magnitute of derivative
%   theta - orientation of derivative (-pi<=theta(i,j)<=pi)
%  Output:
%   edgeMap - a logical map (1 means edge and 0 elsewise)
%
%  Jiani Li, 09/25/2016

[row,col] = size(J);   % size of image

% shifted J map
J_right = zeros(row,col); J_right(:,1:end-1) = J(:,2:end);
J_left = zeros(row,col); J_left(:,2:end) = J(:,1:end-1);
J_up = zeros(row,col); J_up(2:end,:) = J(1:end-1,:);
J_down = zeros(row,col); J_down(1:end-1,:) = J(2:end,:);
J_rightup = zeros(row,col); J_rightup(2:end,1:end-1) = J(1:end-1,2:end);
J_rightdown = zeros(row,col); J_rightdown(1:end-1,1:end-1) = J(2:end,2:end);
J_leftup = zeros(row,col); J_leftup(2:end,2:end) = J(1:end-1,1:end-1);
J_leftdown = zeros(row,col); J_leftdown(1:end-1,2:end) = J(2:end,1:end-1);


us_theta = mod(theta + pi,pi);   % 0<= us_theta < pi, note that the positive direction of theta is clockwise!
% interpolate weight matrix
wy = abs(tan(us_theta)); 
wx = abs(cot(us_theta));

[indx,indy] = meshgrid(1:col,1:row);
not_border = (indx > 1) & (indx < col) & (indy > 1) & (indy < row);

%% non-maximum suppression
% principle when equal: always keep the 'right' and the 'down'

% for pixels: 0 <= us_theta < pi/4 (0 < theta < pi/4 | -pi= < theta < -3pi/4 | theta = pi)
pset = (us_theta >= 0) & (us_theta < pi/4) & not_border;
J1 = (1-wy).*J_right + wy.*J_rightdown;% right & right-up
J2 = (1-wy).*J_left + wy.* J_leftup;% left & left-down
chosen1 = ((J >= J1) & (J >= J2)) & pset;

% for pixels: pi/4 <= us_theta < pi/2 (pi/4 <= theta < pi/2 | -3pi/4 <= theta < -pi/2)
pset = (us_theta >= pi/4) & (us_theta < pi/2) & not_border;
J1 = (1-wx).*J_up + wx.*J_leftup;% right-up & up
J2 = (1-wx).*J_down + wx.*J_rightdown;% left-down & down
chosen2 = ((J >= J1) & (J >= J2)) & pset;

% for pixels: pi/2 <= us_theta < 3pi/4 (pi/2= < theta < 3pi/4 | -pi/2 <= theta < -pi/4 )
pset = (us_theta >= pi/2) & (us_theta < 3*pi/4) & not_border;
J1 = (1-wx).*J_up + wx.*J_rightup;% up & left-up
J2 = (1-wx).*J_down + wx.*J_leftdown;% down & down-right
chosen3 = ((J >= J1) & (J >= J2)) & pset;

% for pixels: 3pi/4 <= us_theta < pi (3pi/4 <= theta < pi | -pi/4 <= theta < 0)
pset = (us_theta >= 3*pi/4) & (us_theta <= pi) & not_border;
J1 = (1-wy).*J_left + wy.*J_leftdown;% left-up & left
J2 = (1-wy).*J_right + wy.*J_rightup;% right-down & right
chosen4 = ((J >= J1) & (J >= J2)) & pset;

chosen = chosen1 | chosen2 | chosen3 | chosen4 ;

edgeMap = (chosen) & (J>0);

end

