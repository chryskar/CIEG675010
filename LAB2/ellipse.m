function [h] = ellipse(ax,xCenter,yCenter,xRadius,yRadius,fig_num)
% function [h] = ellipse(ax,xCenter,yCenter,xRadius,yRadius,fig_num)

%% Ellipse Function
% Author: Chrysostomos Karakasis 702529334

% Description: This function is responsible for plotting an ellipse of
% specific center and x and y radii. The function receives as inputs: the
% axes that the user wishes to be utilized for the plotting of the ellipse,
% the x-coordinate for the center of the ellipse, the y-coordinate for the
% center of the ellipse, the radius of the ellipse in the x-axis, the
% radious of the ellipse in the y-axis and the figure number that the user
% wished to be utilized for the plotting of the ellipse. The function
% returns a variable associated with the plotting of the ellipse for later
% use in a legend function.

theta = 0 : 0.01 : 2*pi; % Angles from 0 to 360 deg

% Polar coordinates will be utilized for the plotting of the ellipse
x = xRadius * cos(theta) + xCenter; % Vector for the x-coordinates of the points on the ellipse
y = yRadius * sin(theta) + yCenter; % Vector for the x-coordinates of the points on the ellipse
figure(fig_num) % Create or open the figure specified from the user
h = plot(ax,x, y,'m','LineWidth', 3); % Plot all points of the ellipse on the axes specified by the user, using the magenta color
axis square; % Make axes square to actually see the ellipse shape
grid on; % Use grid to help user identity distances between objects
end