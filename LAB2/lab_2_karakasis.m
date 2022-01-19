%% CIEG 675 LAB#2 Due Monday January 18, 2021
%% Author: Chrysostomos Karakasis 702529334
close all;
clear all;
% Functions are located in the same folder as this file

%% Part (1)
% Develop a user-defined function that will take as input the mortar firing
% angle, the distance to target and the initial mortar velocity (those 
% better in matlab can also have the target elevation as an input; 
% others:just assume the mortar launch point elevation and target elevation
% are the same at zero). 
angle = 45; % Input mortar firing angle
target_dist = 500; % Distance to target 
v = 70; % Initial mortar velocity
t_elev = 40; % Target elevation
fig_num = 1; % Figure number that the user wishes to be utilized for plotting
% Call mortar function with the aforementioned arguments and receive a
% boolean variable (hit(1) or miss(0)) and the distance achieved between
% mortar and target.
[flag,dist] = mortar(angle, target_dist, v, t_elev,fig_num); 
hold off;
%% Q2 - Circle Problem
% Write a user-defined function(call it circleplotor similar)that will take,
% as input,the coordinates for the center of a circle along withits radius,
% and then plots the circle.It should also plot a marker in the center of 
% the circle (your choice of marker type and color).Make the axis equal so 
% it actually looks like a circle.  As a second step, add,as an input to 
% your function,a color(i.e. ‘k’or ‘r’).  Then use the fill command to fill 
% in the circle with that color.Demonstrate your function works by 
% including 2 plots of 2 different circles (different center coordinates, 
% different radii and different fill colors).

% First plot of a circle
xCenter = 1; % Define the x-coordinate for the center of a circle
yCenter = 2; % Define the y-coordinate for the center of a circle
r = 3; % Define the radius of the circle
col = 'black'; % Define the desired color to be utilized for the filling of the circle
fig_num = 2; % Define the figure that the user wishes to be utilized for the plotting
% Call the circleplot function to plot the circle and receive a "pointer"
% to the plot as a return variable
h2 = circleplot(xCenter,yCenter,r,col,fig_num);

% Second plot of a circle
xCenter = -1; % Define the x-coordinate for the center of a circle
yCenter = -2; % Define the y-coordinate for the center of a circle
r = 6; % Define the radius of the circle
col = 'red'; % Define the desired color to be utilized for the filling of the circle
fig_num = 3; % Define the figure that the user wishes to be utilized for the plotting
% Call the circleplot function to plot the circle and receive a "pointer"
% to the plot as a return variable
h3 = circleplot(xCenter,yCenter,r,col,fig_num);
hold off;

%% Q3 - Make a Plot
% Make a plot of your choosing. Then, using set and get commands, do the
% changes described below.

fig3 = figure(4); % Create a new figure
t = 0:pi/2:100; % Time vector utilized for the plot
r = t.^2/100; % Second variable utilized for the plot
[x,y] = pol2cart(t,r); % Transforms corresponding elements of the polar coordinate arrays theta and rho to two-dimensional Cartesian, or xy, coordinates.
plot(x,y) % Plot a spiral square (why not?)
g1 = get(gcf); % Store the result of the get function so that we don't call it everytime
g2 = get(gca); % Store the result of the get function so that we don't call it everytime

% Increase the figure window width and height to something larger than default
set(gcf,'position',[680 300 700 650])
%Default was [680   558   560   420]

% Change the figure Color to [1 1 1] (white)
set(gcf,'color',[1 1 1])

% Change the XTicks to something other than what they were
set(gca,'XTick',[-95.2 -84 -60.0001 42 55.55 77.77])

% Change the YTicks to something other than what they were
set(gca,'YTick',[-88.2 -74 -50.0001 32 45.55 67.77])

% Change FontSize of tick labels on the x-and y-axis to 14 pt. (i.e. xticklabel)
set(gca,'Fontsize',14)

% Change the FontName to Times New Roman
set(gca,'Fontname','Times New Roman')

% Change the tick direction to out
set(gca,'TickDir','out')

% Add a title, xlabel and ylabel. Make them bold and FontSize 16
set(g2.Title,'string','LAB 2 - Q3 - Figure Title','FontWeight','bold','FontSize',16)
set(g2.XLabel,'string','Figure XLabel','FontWeight','bold','FontSize',16)
set(g2.YLabel,'string','Figure YLabel','FontWeight','bold','FontSize',16)

% Turn the box off (must do this last!)
set(gca,'Box','off')

% Save(print)the figure to a folder other than the current folder WITHOUT
% changing the current folder you are in.  That is, save to a new folder in
% the actual save call using the full filepathname
print(fig3,'C:\Users\chryskar\Desktop\lab_3_Q3','-dpng')
