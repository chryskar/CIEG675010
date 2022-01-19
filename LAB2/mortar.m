function [flag,dist_m_t] = mortar(angle,target_dist,v,t_elev,fig_num)
% function [flag,dist_m_t] = mortar(angle,target_dist,v,t_elev,fig_num)

%% Mortar Function 
% Author: Chrysostomos Karakasis 702529334

% Description: This function is responsible for plotting the trajectory of 
% a fired mortar and determining whether it hits or misses a target. The 
% function will take as inputs the mortar firing angle (in degrees), the 
% distance to target (in meters), the initial mortar velocity (in m/s) and 
% the target elevation (in meters). Also, the function receives as input 
% the figure number that user wishes to be utilized for the plotting of the 
% data. The function will return a boolean variable to inform whether a hit
% (1) or a miss (0) occured. In the latter case, the minimum 2D-distance 
% achieved between the target and the mortar will be also returned, while 
% in the former case the distance between the mortar and the target at the 
% time of impact is returned. We assume that the target elevation is 
% stricly non-negative and that if the mortar enters a 10 meter radius away
% from the target, then we have a hit. In the case where the mortar does 
% not enter the 10 meter radius, then it will land when the y-coordinate 
% becomes zero. Moreover, we neglect air resistance, Coriolis, earth 
% surface curvature, etc.

% We will define two sets of axes. The first one is for plotting the
% trajectory of the mortar 
ax1 = axes('Position',[0.1 0.1 0.7 0.7]);  

angle = deg2rad(angle); % Change the mortar firing angle to radians
g = 9.81; %m/s - Define gravitational acceleration

figure(fig_num) % Create a new figure based on the number given by the user

%% Draw Ellipse for target
xCenter = target_dist; % The x-coordinate of the target is equal to the distance of the target from the launch location
yCenter = 2+t_elev; % The y-coordinate of the target is equal to the target elevation plus a small offset to have the target lying completely overground
xRadius = 5; % The radius of the ellipse in the x-axis around the target 
yRadius = 2; % The radius of the ellipse in the y-axis around the target
h1 = ellipse(ax1,xCenter,yCenter,xRadius,yRadius,fig_num); % Use custom-defined function to draw ellipse around target

%% Draw 10 meter radius around the target
h2 = viscircles(ax1,[xCenter yCenter],10,'Color','r','LineStyle','--'); % We use this function to create a circle around the target at 10m
title('Mortar Trajectory Problem - Q1','interpreter','latex','FontSize',14); % Add a temporary title

%% Plotting of the trajectory
x_coord = 0; %m Position of the launch site in the x-axis
y_coord = 0; %m Position of the launch site in the y-axis

t = 0; % Counter of iterations in the while loop utilized as time vector
flag = 0; % Boolean variable that will specify whether a hit (1) or a miss (0) occured
x_coord_v = []; % Vector where the x-coordinates of all positions of the mortar throughout the trajectory will be saved
y_coord_v = []; % Vector where the y-coordinates of all positions of the mortar throughout the trajectory will be saved
min_dist = Inf; % Initialization of the minimum distance achieved between the mortar and the target

while y_coord >= 0 % As long as the mortar has a nonnegative height, the simulation goes on
    figure(fig_num); % Keep plotting on the figure specified by the user
    axis(ax1,'equal'); % Make axis equal so that we observe the ellipse and circle shapes
    hold on; % Retains plots in the current axes so that new plots added to the axes do not delete existing plots
    h3 = plot(ax1,x_coord,y_coord,'xb','LineWidth',2); % Plot currect position of mortar using the 'x' marker and the 'blue' color
%     title('Mortar Trajectory Problem - Q1','interpreter','latex','FontSize',14); % Add a temporary title
    pause(0.1); % Pause for 0.1 seconds to see the intermediate positions
    
    t = t+0.1; % Increase iteration-time counter
    x_coord = v*cos(angle)*t; % Calculate new x-coordinate of the mortar's updated position in meters
    y_coord = v*sin(angle)*t - 0.5*g*t.^2; % Calculate new y-coordinate of the mortar's updated position in meters
    x_coord_v = [x_coord_v x_coord]; % Store new x-coordinate to the vector with previous values
    y_coord_v = [y_coord_v y_coord]; % Store new y-coordinate to the vector with previous values   
    
    distance_v = [xCenter yCenter; x_coord y_coord]; % Current location of mortar and target
    dist_m_t = pdist(distance_v); %2D-Distance(euclidean) between the target and the mortar
    if dist_m_t < min_dist % Compare with minimum distance found so far
        min_dist = dist_m_t; % If current distance is less than the so far minimum, then update minimum value
    end
    
    if dist_m_t <= 10 % If current distance is less than 10 meters, then we have a hit
        flag = 1; % Update boolean variable in case of a hit
        break; % Since a hit occured there is not point of continuing the simulation
    end
end

legend([h1 h2 h3(1)],'Target','10m Radius','Trajectory Points of Mortar','interpreter','latex','FontSize',14); % Add labels for everything plotted in the figure
xlabel('Horizontal Distance (m)','interpreter','latex','FontSize',14); % Add a label for the x-axis
ylabel('Vertical Distance (m)','interpreter','latex','FontSize',14); % Add a label for the y-axis

fig1 = figure(fig_num); % Get a variable associated with fig_num to refer to it
set(gcf, 'Position', get(0, 'Screensize')); % Change the figure to full-screen
if flag == 1 % In case of a HIT
    disp(sprintf(['HIT! Target was destroyed!\nDistance at impact: ',num2str(min_dist),' (m)'])); % Display message to the command window to inform the user about the hit event and the mortar-target distance at impact
    text(ax1,mean(fig1.CurrentAxes.XLim),mean(fig1.CurrentAxes.YLim),'HIT!','Color','red','FontSize',20,'HorizontalAlignment','center'); % Add text to the center of the first axes to inform the user about the result of the shot
    title({['Mortar Trajectory Problem - Q1'],['Distance between mortar and target at impact: ',num2str(min_dist),' (m)']},'interpreter','latex','FontSize',14); % Update the title of the figure with the mortar-target distance at impact
else
    disp(sprintf(['MISS! You missed the target...\nMinimum Distance: ',num2str(min_dist),' (m)\nTry again!'])); % Display message to the command window to inform the user about the miss event and the minimum mortar-target distance achieved
    text(ax1,mean(fig1.CurrentAxes.XLim),mean(fig1.CurrentAxes.YLim),'MISS!','Color','red','FontSize',20,'HorizontalAlignment','center'); % Add text to the center of the first axes to inform the user about the result of the shot
    title({['Mortar Trajectory Problem - Q1'],['Minimum Distance achieved between mortar and target: ',num2str(min_dist),' (m)']},'interpreter','latex','FontSize',14); % Update the title of the figure with the minimum distance achieved
end

%% Define second set of axes to zoom over the target
ax2 = axes('Position',[0.65 0.2 0.28 0.28],'Box','on'); % Define second set of axes
ellipse(ax2,xCenter,yCenter,xRadius,yRadius,fig_num); % Draw again the ellipse around the target at the new axes
hold on; % Retains plots in the current axes so that new plots added to the axes do not delete existing plots
viscircles(ax2,[xCenter yCenter],10,'Color','r','LineStyle','--'); % Draw again the 10 meter circle around the target at the new axes
plot(ax2,x_coord_v,y_coord_v,'xb','LineWidth',2); % Plot all positions of the mortar trajectory
axis(ax2,[xCenter-15 xCenter+15 yCenter-15 yCenter+15]); % Define the axis limits of the second set of axes to observe a close neighborhood of the target
end