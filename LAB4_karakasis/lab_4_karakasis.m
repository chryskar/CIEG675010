%% CIEG 675 LAB #4 Due Tuesday February 1, 2021
%% Author: Chrysostomos Karakasis 702529334
close all;
clear all;
% Functions are located in the same folder as this file

%% Part (1)
% The goal at the end of this problem is to create a .avimovie file that is
% a time lapse of the beginning of the Blizzard in Newark, DE of January 2016.

% Download the compressed folder from Canvas called "BlizzardImages.zip" and
% extract the images to your computer. You will need to write a script to
% find all of the .JPG files in the directory where you extracted the
% images. Recall the dir command discussed in class.
D = dir('./BlizzardImages/'); % We assume that a folder exists inside the folder that contains this code

% Then, with a for loop, write a code that will create a time lapse of the
% blizzard (using imshow) and write a movie to a .avifile. You should also
% add text somewhere on each image in the time lapse (use a large fontsize)
% that displays the total amount of time that has elapsed, in minutes(you may
% have to refresh yourself on dates in MATLAB and converting dates from
% string format to the number of days relative to 0 AD).In addition, you
% may want to try playing around with the figure window size/position, as
% well as the position of the image inside the figure window, if you want
% to get the movie to look more professional without white/grey borders
% around the image.

% Define the name of the .avi video
vidObj = VideoWriter('Karakasis_Chrysostomos.avi');
% Use a framerate of 8 frames per second (FPS).
FrameRate = 8;
% Open the .avi file
open(vidObj);

for i = 3:length(D) % Skip the two first lines ('.' and '..')
    folder = D(i).folder; % Location of the folder containing the image
    image = D(i).name; % Name of image
    path = [folder '\' image]; % Location of the image
    imshow([D(i).folder '\' image]); % Show the image
    set(gcf,'Position',[446   233   933   700]); % Change the size of the image to minimize any gray areas
    
    % maximize the image within a figure window
    im = get(gcf,'children');
    set(im,'units','normalized','position',[0 0 1 1]);
    
    % Use the 'day2secthat' function to find elapsed time in seconds and
    % divide by 60 to get elapsed time in minutes
    time_elap = day2secthat([D(3).datenum D(i).datenum])/60;
    
    % Print time elapsed on image
    dim = [0.04 .675 .3 .3]; % Location of the text box in the figure
    str = ['+',num2str(time_elap(2),'%6.2f') ' minutes elapsed']; % The text that will appear on the figure
    annotation('textbox',dim,'String',str,'FitBoxToText','on','FontSize',16,'color','y','EdgeColor','y');
    
    pause(0.01) % so we can see the change while making it
    currframe = getframe(gcf); % "grab" the current frame
    writeVideo(vidObj,currframe); % write the current frame
end
close(vidObj); % Close the .avi file

%% Part (2)
% Write a function (call it gridinterp) to do an inverse distance weighting
% algorithm for 2D interpolation. It is a technique for interpolating to a
% uniform grid.  Your function should have 7 input variables: x,y, z, xg,
% yg, d and alpha. Your function should have a single output vector, zg,

% Load the "Avalon_survey.mat" file that we will use for testing
avalon = load('Avalon_survey.mat');
% Save the scatted x data
x = avalon.x;
% Save the scatted y data
y = avalon.y;
% Save the scatted z data
z = avalon.z;

% Load default values for the three parameters
dx = 4; % Define increment for the desired uniform grid
d = 10; % Define the search radius
alpha = 1; % 1 <= a <= 2 order to inversely weight each point at x,y

xx=-10:dx:100; % x locations of the desired uniform grid
yy=-180:dx:110; % y locations of the desired uniform grid
[X,Y]=meshgrid(xx,yy); % vectorize the desired uniform grid
xg=X(:); % horizontal grid location in x-axis
yg=Y(:); % horizontal grid location in y-axis
% Basically xx is repeated 'length(yy)' times and yy is repeated
% 'length(xx) times, to create a grid

zg=gridinterp(x,y,z,xg,yg,d,alpha); % Call the function we created to apply the interpolation
ZG=reshape(zg,size(X)); % Reshape zg from a vector to 2D-matrix with appropriate size to use surf plot

fig2 = figure(2); % Create a new figure
clf; % Clear anything that might exist already in the figure
surf(X,Y,ZG); % 3D-plot the interpolated data wrt to the desired uniform grid
colormap(jet); % Change to a RGB spectrum
colorbar % Add a colorbar to have an indication of the values at the surface's points
xlabel('cross-shore distance (m)','fontsize',14,'interpreter','latex'); % Set an appropriate x-axis
ylabel('along-shore distance (m)','fontsize',14,'interpreter','latex'); % Set an appropriate y-axis
zlabel('Elevation (m)','fontsize',14,'interpreter','latex'); % Set an appropriate z-axis
% Set a title to show which values were utilized for (dx,d,alpha)
title({'LAB \# 4 - Part (2) - gridinterp',['dx = ',num2str(dx),', d = ',num2str(d),', $\alpha$ = ',num2str(alpha)]},'fontsize',14,'interpreter','latex');
set(gcf, 'Position', get(0, 'Screensize'));  % Switch to full-screen to have a better visualization of the figure
az = 47;
el = 39;
view(az,el); % Change the camera view angle to match the one in the given example
print(fig2,'lab_4_part_2_default','-depsc','-r600'); % Save figure in a colored .eps format using 600 dpi resolution

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Explore effect of parameter dx, d, alpha
% Defaults: dx=4, d=10, alpha=1
dx_v = [4 10 20];
d_v = [1 10 100];
alpha_v = [1 1.5 2];
for i = 1:3 % Repeat for dx, d and alpha
    for j = 1:3 % Try three different values for each parameter
        if i == 1 % Change only dx paremeter and load default for d and alpha
            param = 'dx'; % Define which parameter is currently explored
            dx = dx_v(j); % Define increment for the desired uniform grid
            d = 10; % (Default value) Define the search radius
            alpha = 1; % (Default value) 1 <= a <= 2 order to inversely weight each point at x,y
        elseif i == 2 % Change only d paremeter and load default for dx and alpha
            param = 'd'; % Define which parameter is currently explored
            dx = 4; % (Default value) Define increment for the desired uniform grid
            d = d_v(j); % Define the search radius
            alpha = 1; % (Default value) 1 <= a <= 2 order to inversely weight each point at x,y
        else % Change only alpha paremeter and load default for dx and d
            param = 'alpha'; % Define which parameter is currently explored
            dx = 4; % (Default value) Define increment for the desired uniform grid
            d = 10; % (Default value) Define the search radius
            alpha = alpha_v(j); % 1 <= a <= 2 order to inversely weight each point at x,y
        end
        
        xx=-10:dx:100; % x locations of the desired uniform grid
        yy=-180:dx:110; % y locations of the desired uniform grid
        [X,Y]=meshgrid(xx,yy); % vectorize the desired uniform grid
        xg=X(:); % horizontal grid location in x-axis
        yg=Y(:); % horizontal grid location in y-axis
        % Basically xx is repeated 'length(yy)' times and yy is repeated
        % 'length(xx) times, to create a grid
        
        zg=gridinterp(x,y,z,xg,yg,d,alpha); % Call the function we created to apply the interpolation
        ZG=reshape(zg,size(X)); % Reshape zg from a vector to 2D-matrix with appropriate size to use surf plot
        
        fig3(i) = figure(2+i); % Create a new figure (one for each variable)
        hold on;
        sgtitle({['LAB \# 4 - Part (2) - gridinterp','Explore the parameter ``',param,'"']},'fontsize',14,'interpreter','latex');
        s = subplot(3,1,j);
        surf(X,Y,ZG) % 3D-plot the interpolated data wrt to the desired uniform grid
        title(['dx = ',num2str(dx),', d = ',num2str(d),', $\alpha$ = ',num2str(alpha)],'fontsize',14,'interpreter','latex');
        az = 47;
        el = 39;
        % Change the camera view angle to match the one in the given example
        view(s,az,el);
        
        colormap(jet); % Change to a RGB spectrum
        colorbar; % Add a colorbar to have an indication of the values at the surface's points
        xlabel('cross-shore distance (m)','fontsize',14,'interpreter','latex'); % Set an appropriate x-axis
        ylabel('along-shore distance (m)','fontsize',14,'interpreter','latex'); % Set an appropriate y-axis
        zlabel('Elevation (m)','fontsize',14,'interpreter','latex'); % Set an appropriate z-axis
        set(gcf, 'Position', get(0, 'Screensize'));  % Switch to full-screen to have a better visualization of the figure
    end
    % For each variable save the corresponding figure with a relevant name
    print(fig3(i),['lab_4_part_2_' param 's'],'-depsc','-r600'); % Save figure in a colored .eps format using 600 dpi resolution
end

%% Part (3)
% OPTIONAL: COVID-19 can be modeled with the SIR model. The model consists
% of three coupled differential equations
for b1 = 1 : 10 % Try 10 different values
    k1 = 14; % Define the average time to recover
    
    S_0 = 1; % Initial susceptible rate
    I_0 = 0.01; % Initial infected rate (very small)
    R_0 = 0; % Initial recovered rate
    figure_num = 6;
    max_days = 200; % Solve the model for 200 days
    [p1,p2,p3] = covid_sir(S_0,I_0,R_0,k1,b1,max_days,figure_num);
end

% Set a colorbar to indicate the value of "b1" at each curve
% Yellow corresponds to b1=1 and red corresponds to b1=10
h = colorbar;
set(h,'Limits',[0.1 1]); % Keep colors from yellow to red
colormap([ones(10,1) [1:-0.1:0.1]' zeros(10,1)]);
set(h,'ticklabels',[1:1:10]); % Change the tick labels to show the values of parameter b1

% Set a corresponding legend
legend([p1(end) p2(end) p3(end)],'Susceptible (S)','Infected: (I)','Recovered/Removed: (R)','interpreter','latex');
xlabel('Time (Days)','interpreter','latex'); % x-axis (time in days)
ylabel('Normalized Population','interpreter','latex'); % y-label (normalized population)
ylabel(h, 'Values for parameter ``$b_{1}$"','interpreter','latex','Fontsize',14); % ylabel of colorbar (values of b1)
ylim([0 1.3]) % Increase the size to avoid overlapping with legend
% Set appropriate title to explain this figure
title({'LAB \# 4 - Part (3) - COVID-19 SIR Model','Explore the parameter ``$b_{1}$"'},'fontsize',14,'interpreter','latex');
print(figure(figure_num),'lab_4_part_3','-depsc','-r600'); % Save figure in a colored .eps format using 600 dpi resolution

%% Part (4)
% Solve dy/dx = -y(x) subject to initial condition y(0)=1. Note the real
% answer is exp(-x). Do the solution in matlab two different ways.
% A) Solve using a forward difference.
% B) Solve using a central difference.
% For each solution approach, use three different dx values and make a
% table showing the root-mean-square error (compared to the real answer) as
% a function of dx and numerical scheme. Also, show a plot of your
% solutions.

% OPTIONAL: Solve the same equation by exploring the use of the functions
% ODE23 and/or ODE45.  These are matlab's built in differential equation
% solvers (there are others) for differential equations over a certain
% range in independent variable and with a known boundary condition.  Make
% a plot of your solutions from 1a (A and B) with the solution from one of
% these ODE solvers.
clear all; % Delete any previous data

% Part (4.a) - Forward Difference
% df(x)/dx = {f(x+Dx)-f(x)}/Dx => y(x+Dx) = y_dot(x)*Dx + y(x)

% Part (4.b) - Central Difference
% df(x)/dx = (1/2)*{f(x+Dx)-f(x-Dx)}/Dx

% In this case you wish to solve this equation y_dot(x) = -y(x)
dx_v = [0.5 0.1 0.001]; % Vector containing the three values for dx

fig7 = figure(7); % Open a new figure
clf; % Delete everything in this figure in case it was open before execution

for i = 1:3 % Repeat the same process for all three different dx values
    clear x y_f y_c; % Clear every time the variables that will be utilized later
    dx = dx_v(i); % Load the desired dx value for this iteration (Step size)
    endx = 4; % Where you want to stop calculations
    y_f(1) = 1;  % the initial value in the y vector corresponding to x=0 (forward case)
    y_c(1) = 1;  % the initial value in the y vector corresponding to x=0 (central case)
    
    x(1) = 0; % initial x location
    ct = 1;  % initialize counter
    
    while x(ct)<endx % Keep solving the ODE until reaching the endx threshold
        ct = ct+1; % increment counter
        x(ct) = x(ct-1)+dx;   % increment x (keep it to simulate real solution later)
        
        % Forward Difference: df(x)/dx = {f(x+Dx)-f(x)}/Dx => y(x+Dx) = y_dot(x)*Dx + y(x)
        y_f(ct) = y_f(ct-1)-dx*y_f(ct-1);  % find the new y for the new x using the forward difference
        
        % Central Difference: df(x)/dx = (1/2)*{f(x+Dx)-f(x-Dx)}/Dx =>
        % y_dot(x)*Dx*2 = y(x+Dx)-y(x-Dx)=> y(x+Dx) = y(x-Dx) + y_dot(x)*Dx*2
        y_c(ct) = y_c(ct-1)+2*dx*(-y_c(ct-1)); % find the new y for the new x using the central difference
    end
    yreal = exp(-x);  % this is the real solution found from integration
    
    subplot(3,1,i); % Plot the results for all three cases
    hold on; % Plot everything in the same subplot
    plot(x,yreal,'k'); % Plot the real solution
    plot(x,y_f,'ro'); % Plot the approximate solution for the forward case
    plot(x,y_c,'bd'); % Plot the approximate solution for the central case
    legend('real soln','approximation - forward','approximation - central','interpreter','latex'); % Set appropriate legend
    title({['$dx=$',num2str(dx)],['Forward RMS: ',num2str(rms(yreal-y_f)),' - Central RMS: ',num2str(rms(yreal-y_c))]},'interpreter','latex');
    sgtitle({'LAB \# 4 - Part (4) - Numerical Methods','Explore the parameter ``$dx$"'},'fontsize',14,'interpreter','latex');
    xlabel('x');
    ylabel('y');
end

% Time for the optional part

% ODE Solver
syms y(t) [1 1] real; % Define one symbolic variable
% Plug in the desired ODE
eqn1 = diff(y1(t),t) == -y1; % Define the equation we wish to solve

eqns = [eqn1]; % (Optinal step, mostly useful when we wish to solve numerous coupled equations like in Part (3))
vars = [y1(t)]; % (Optinal step, mostly useful when we wish to define numerous coupled variables like in Part (3))
[M,F] = massMatrixForm(eqns,vars); % Retrieve the mass matrix form (might be an overkill but I am used to this methodology)
f = M\F; % Retrieve the system that we wish to solve in a compact form

% Set the desired initial conditions for the ODE
init_conds = [1];

odefun = odeFunction(f,vars); % Plug in everything using the odefun function to utilize afterwards the ODE45 solver
t_total=[]; % In this vector we will store all time variables that the ode45 solver returns
sol=[];     %In this vector we will store all state variables that the ode45 solver returns

i=1;        % Iteration counter for the outer loop
i_th = 100; % Maximum Number of Iterations
opts = odeset('RelTol',1e-6,'AbsTol',1e-6); %Absolute and Relative error tolerances for ode45
max_step_size = dx_v(3); % Utilize the best combination from above
while(1) % Keep solving the ODE until reaching the desired number of time steps
    [t_i,sol_i] = ode45(odefun,[0 max_step_size],init_conds,opts); %Solve the ODE using the specified ICs for dx time steps with the error tolerances "opts" specified earlier
    t_total=[t_total;t_i+((i-1)*max_step_size)]; % Adjust the local time vector, so that it corresponds to the total time elapsed and append it to the total time vector
    sol=[sol;sol_i]; % Append the local solutions for the state to the total state solutions vector
    init_conds=[sol(end,1)]; % Set the last state of the system for this iteration, as ICs for the next one
    i = i + 1; % Increase the iteration counter
    if t_total(end) > endx % Check whether the system has reached the desired number of time steps
        break % Break in case the desired number of time steps has been reached
    end
end

yreal_ode = exp(-t_total);  % this is the real solution found from integration (for this time step)

fig8 = figure(8); % Open a new figure
clf; % Clear everything in this figure
hold on; % Plot everything in the same figure
plot(x,yreal,'k'); % Plot the real solution
plot(x,y_f,'ro'); % Plot the approximate solution using the forward difference (was more accurate than central)
plot(t_total,sol,'bd'); % Plot the ode45 solution
legend('real soln','approximation - forward','ode45','interpreter','latex'); % Set appropriate legend
title({'Compare Forward Difference to ODE45',['Forward RMS: ',num2str(rms(yreal-y_f)),' - ODE45 RMS: ',num2str(rms(yreal_ode-sol))]},'interpreter','latex');
xlabel('x');
ylabel('y');

set(fig7,'units','normalized','Position',[0.2641    0.3657    0.3812    0.5389]); % Increase the figure size a bit.
print(fig7,'lab_4_part_4_dxs','-depsc','-r600'); % Save figure in a colored .eps format using 600 dpi resolution
print(fig8,'lab_4_part_4_for_vs_ode45','-depsc','-r600'); % Save figure in a colored .eps format using 600 dpi resolution

%% Part (5)
% Load the image speed_limit.JPG from canvas and explore the matlab
% function called 'edge' to perform edge detection on the simple image.
% Obtain the output from edge and change the background image color to
% yellow and the letter and number edges to blue. There are multiple ways
% to do this and you will need to figure out how to go from the original
% RGB image to grayscale for analysis and then back to RGB for output.

% Load the image which should located in the same folder as this code
im = imread('speed_limit.JPG');

% Convert image from RGB to grayscale
I = rgb2gray(im);

% Detect edges using the Sobel method
BW1 = edge(I,'sobel');

% After the edge detection, matrix BW1 will contain for each element of the
% image either a 0 or a 1, in order to indicate whether that that
% element-pixel corresponds to an edge or not (background)

% Turn all zeros to yellow( = [1 1 0] in RGB) - background elements
[back_row,back_col] = find(BW1 == 0);
% Turn all ones to blue(=[0 0 1] in RGB) - letter/numbers edges
[let_num_row,let_num_col] = find(BW1 == 1);

% In order to change the colors of the image in the RGB format, we have to
% specify for each element the corresponding RGB combination that
% corresponds to the desired color. For instance, im(1,1,1) corresponds to
% the R(red) percentage of the first "pixel", im(1,1,2) corresponds to the
% G(green) percentage of the first "pixel" and finally im(1,1,3)
% corresponds to the B(blue) percentage of the first "pixel". Hence, the
% above indices correspond to the location of every pixel that we will to
% change its color.

for i = 1:length(back_row) % Go through all background elements-pixels
    im(back_row(i),back_col(i),1) = 255; % Change the R(red) part to 1(255)
    im(back_row(i),back_col(i),2) = 255; % Change the G(green) part to 1(255)
    im(back_row(i),back_col(i),3) = 0; % Change the B(blue) part to 0(0)
    % Hence, now all background pixels have the color [1 1 0] which corresponds
    % to yellow
end
for i = 1:length(let_num_row) % Go through all edge elements-pixels
    im(let_num_row(i),let_num_col(i),1) = 0; % Change the R(red) part to 0(0)
    im(let_num_row(i),let_num_col(i),2) = 0; % Change the G(green) part to 0(0)
    im(let_num_row(i),let_num_col(i),3) = 255; % Change the B(blue) part to 255(1)
    % Hence, now all edge pixels have the color [0 0 1] which corresponds
    % to blue
end

figure(9); % Open a new figure
imshow(im); % Show the modified image
imwrite(im,'speed_limit_mod.JPG'); % Save the modified image