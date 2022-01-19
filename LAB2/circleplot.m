function [h] = circleplot(xCenter,yCenter,r,col,fig_num)
% function [h] = circleplot(xCenter,yCenter,r,col,fig_num)

%% Circleplot Function
% Author: Chrysostomos Karakasis 702529334

% Description: This function is responsible for plotting a circle of
% specific center and radius. Furthermore, it receives as an input a color
% (i.e. ‘k’or ‘r’), which is utilized to fill in the circle with that color.
% It also plots a marker in the center of the circle. Finally, the function
% returns a variable linked to the generated plot.

angle_v = 0:0.01:2*pi; %Angles from 0 to 360 deg

%We will use polar coordinates to create the circles
x1 = xCenter + r*cos(angle_v); %X-coordinates of the points of the first circle derived from the polar coordinates 
y1 = yCenter + r*sin(angle_v); %Y-coordinates of the points of the first circle derived from the polar coordinates

figure(fig_num) % Create or open the figure specified by the user
axis equal % Make the axis equal so it actually looks like a circle.
hold on; %Use "hold on" to plot all three circles on the same plot

fill(x1,y1,col,'linewidth',3) %Plot the circle first, fill in the inside with the desired color (red) and set line width to 3
h = plot(xCenter,yCenter,'b+','markersize',10); % Plot the center of the circle and add a "blue" "+" marker 

% Add a title for the plot containing information about the figure.
title({['Example of ``circleplot" Function'],['Center Coordinates: (',num2str(xCenter),', ',num2str(yCenter),') - Radius: ',num2str(r),' - Fill Color: ',col]},'interpreter','latex')
legend('Circle filled with color','Center of Circle','interpreter','latex') % Add a legend specifying everything in the figure
end