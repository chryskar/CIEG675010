function zg = gridinterp(x, y, z, xg, yg, d, alpha)
% function zg = gridinterp(x, y, z, xg, yg, d, alpha)

% Description:
% This function performs an inverse distance weighting algorithm for 2D 
% interpolation. It is a technique for interpolating to a uniform grid. It
% receives as inputs a measured parameter ('z') and the horizontal locations 
% ('x') and ('y') at which it is scattered at. Furthermore, it receives the
% horizontal locations ('xg' and 'yd') of the desired uniform grid that the
% user wishes to utilize for the interpolation. Besides that, the function
% receives as an input a parameter ('d'), which specifies the allowable
% radial distance from each uniform grid point that the algorithm should
% search for measured data point over which to interpolate (search radius).
% Finally, the last input of the function is parameter a(alpha), which is 
% the order to inversely weight each point at x, y based on its distance 
% from xg, yg (i.e., for alpha = 2, the weighting is inverse distance 
% squared). The output of the function is an interpolated version of 'z' 
% onto the desired uniform grid given by horizontal grid locations xg and
% yg.

% INPUTS:
%x: horizontal location of measured parameter
%y: vertical location of measured parameter
%z: measured parameter (e.g. elevation)
%xg: horizontal location of desired uniform grid
%yg: vertical location of desired uniform grid
%d: allowable radial distance from each uniform grid point(search radius)
%alpha: order to inversely weight each point at x, y based on its distance from xg, yg(i.e., for alpha= 2, the weighting is inverse distance squared). Hint, xgand ygshould come from "vectorizing" the uniform grid you develop using the meshgrid function (see below).

% OUTPUTS:
% zg: interpolated version of input variable 'z' onto the desired uniform grid given by xg and yg 

zg = []; % Initialization of the output vector
for j = 1 : length(xg) % Iterate for all desired points
    sum_num = 0; % Initialize the sum that corresponds to the numerator of the final function
    sum_den = 0; % Initialize the sum that corresponds to the denominator of the final function
    for i = 1 : length(x)  % Iterate through all initial points 
        weight_num = norm([xg(j)-x(i) yg(j)-y(i)]); % Weights - Euclidean distance from the grid point to the real value
        if weight_num < d % Apply the search radius and determine whether a point should be taken into account or not
            sum_num = sum_num + z(i)/(weight_num^alpha); % Update the numerator sum in case a point should be included
            sum_den = sum_den + 1/(weight_num^alpha); % Update the denominator sum in case a point should be included
        end
    end
    zgj = sum_num/sum_den; % Apply the final formula and calculate another element of the interpolated output
    zg = [zg; zgj]; % Append the new element to the previous elements of the interpolated output
end

end