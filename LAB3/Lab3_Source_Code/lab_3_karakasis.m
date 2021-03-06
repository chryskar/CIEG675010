%% CIEG 675 LAB#3 Due Tuesday January 26, 2021
%% Author: Chrysostomos Karakasis 702529334
close all;
clear all;
% Functions are located in the same folder as this file

%% Part (1)
% Write a function that will take an arbitrary, but smooth, vector of data
% and locate ALL the local maxima and minima in the vector (call it min_max,
% as "minmax" is a built-in function in newer MATLAB versions starting with
% R2014b). It should return, as output, the indices of the input vector
% where these maxima and minima occur. Test your function on the following
% data set that you should make quite smooth by using small increments in
% x, where x should go from 0 to 10:
% y = x^{1.01}+4*cos(3*pi*x/4)-2*sin(2*pi*x/3)-0.25
% In the same figure, plot y vs x, with the maxima and minima (computed by
% your function, min_max) on top of the curve as different symbols
% (remember you can use >> help plotto see a list of available line/marker
% specifiers). DO NOT USE a built in function to find local minima or
% maxima. I want you to write the code to do it.

x = 0:0.01:10; % Create the x-vector with a small enough increment
y = x.^(1.01)+4*cos(3*pi*x/4)-2*sin(2*pi*x/3)-0.25; % Calculate the smooth dataset that will be used as a testing input for our function
[min_ind,max_ind] = min_max(y) % Call the custom-made 'min_max' to get the locations of the local extrema using the first-derivative test
fig1 = figure(1); % Open a new figure
hold on; % Use "hold on" to plot all three parabolas on the same plot
plot(x,y,'LineWidth',2) % Plot the y vs x data
plot(x(min_ind),y(min_ind),'ro','LineWidth',2) % Plot the local minima using a red 'o' symbol
plot(x(max_ind),y(max_ind),'mx','LineWidth',2) % Plot the local maxima using a magenta 'x' symbol
% Create a new legend specifying everything on the figure
legend('Arbitrary-smooth Vector of Data','Local Minima','Local Maxima','interpreter','latex','location','southeast','FontSize',11)
title('LAB \#3 - Part (1) - Min\_max','FontSize',14,'interpreter','latex') % Create a title for the figure
xlabel('X-axis','FontSize',14,'interpreter','latex') % Set a label for the x-axis
ylabel('Y-axis','FontSize',14,'interpreter','latex') % Set a label for the y-axis
axis([0 10 -10 14]) % Increase axis limits
print(fig1,'lab_3_part_1','-depsc','-r600'); % Save figure in a colored .eps format using 600 dpi resolution

%% Part (2)
% Create a 1 x 2 structure array with 3 fields of your choosing. The only
% stipulation is that the 1st field must be a string (class: char), the
% 2nd field must be a cell array (class: cell), and the 3rd field should be
% an array of numeric values (class: double).

D(1,1).field1 = 'hey'; % Define the first field of the first element of the struct structure (char:string)
D(1,1).field2 = {2}; % Define the second field of the first element of the struct structure (cell)
D(1,1).field3 = [3.24:0.1:5.46]; % Define the third field of the first element of the struct structure (array of doubles)

D(1,2).field1 = 'there'; % Define the first field of the second element of the struct structure (char:string)
D(1,2).field2 = {4}; % Define the second field of the second element of the struct structure (cell)
D(1,2).field3 = [5.46:-0.1:3.24]; % Define the third field of the second element of the struct structure (array of doubles)

D % Show that D is a 1x2 structure array
class(D(1).field1) % Show that the 1st field is indeed a string (class: char)
class(D(1).field2) % Show that the 2st field is indeed a cell array (class: cell)
class(D(1).field3) % Show that the 1st field is indeed an array of doubles (class: double)

%% Part (3)
% Create a 12 x 1 cell array whose contents are the twelve months of the year, in order, starting with 'January'.
c_a_3 = {'January','February','March','April','May','June','July','August',...
    'September','October','November','December'}

%% Part (4)
% Write a function called day2secthat takes an input vector of MATLAB dates
% and converts it first to time relative to the 1st date value in the input
% vector (i.e. the first value in the output vector should be zero) and
% then from days to seconds. The output vector should be the number of
% seconds that have elapsed since the 1st value in the input vector.

% First, we will create an array with four arbitrary dates in an ascending
% order and then we will convert them to MATLAB dates using the 'datenum'
% built-in function
mat_dates_v = [datenum(1973,1,30,15,32,11) datenum(1973,1,30,15,32,31) datenum(1995,10,6,12,30,24) datenum(2021,1,21,13,18,46)];
num_sec = day2secthat(mat_dates_v) % Call the 'day2secthat' function

% Two versions of the function 'day2secthat' were created
num_sec_2 = day2secthat2(mat_dates_v) % Call the 'day2secthat2' function

%% Part (5)
% Download the file from Canvas called "AtlanticCity_TemperatureData.csv".
% Load the data into MATLAB using dlmread(or some other function). The 1st
% column contains dates in MATLAB time. The 2nd column is air temperature
% and the 3rd column is water temperature, both in degrees Celsius. The data
% span all of 2015. Subtract the first date value from all other dates,
% then plot air temperature and water temperature versus time, on the same
% axes. Add a legend, labels, etc. Finally, change the 'xticklabel' values
% of the current axes to be your cell array from Problem 3, where the
% 'xtick' locations should be the number of days after January 1st, for each
% first day of the month (i.e. [1, 32, 60, 91, 121, 152, 182, 213, 244, 274,
% 305, 335]). Resize the figure to be wider, so that all the month labels
% on the x-axis are not overlapping.

% Load the data into MATLAB using dlmread
atl_temp_data = dlmread('AtlanticCity_TemperatureData.csv');

% Create a copy of the loaded data
atl_temp_data_norm = atl_temp_data;
% Subtract the first date value from all other dates
atl_temp_data_norm(:,1) = atl_temp_data_norm(:,1) - atl_temp_data_norm(1,1);

% Plot air temperature and water temperature versus time, on the same axes
fig2 = figure(2);
hold on; % Use "hold on" to plot all three parabolas on the same plot
plot(atl_temp_data_norm(:,1),atl_temp_data_norm(:,2),'LineWidth',1); % Plot the air temparature data wrt time
plot(atl_temp_data_norm(:,1),atl_temp_data_norm(:,3),'LineWidth',1); % Plot the water temparature data wrt time
% Create a new legend specifying everything on the figure
legend('Air Temperature ($C^{\circ}$)','Water Temperature ($C^{\circ}$)','interpreter','latex','FontSize',14,'location','southeast')
title('LAB \#3 - Part (5) - Atlantic City Temperature Data','FontSize',14,'interpreter','latex') % Create a title for the figure
xlabel('X-axis','FontSize',14,'interpreter','latex') % Set a label for the x-axis
ylabel('Y-axis','FontSize',14,'interpreter','latex') % Set a label for the y-axis

% In order to isolate each first day of the month. We assume that for every
% first day of the month, the first entry data will be utilized.
temp = datestr(atl_temp_data(:,1)); % Convert the MATLAB dates back to string format
first_days_per_month = []; % Initialize the array that will store the indices depicting the first entry for the first day of every month

for i = 1:size(atl_temp_data,1) % Iterate through all dates
    if temp(i,1:2)=='01' % Check whether an element corresponds to the first day of the month
        if (temp(i,end-7:end)=='00:00:00') % Check whether an element also corresponds to the first entry of the day
            first_days_per_month = [first_days_per_month i]; % When a first-entry/first-day-month is found, append it to the existing vector
        end
    end
end

% Change the 'xticks' locations to the number of days after January 1st, for each
% first day of the month.
xticks(atl_temp_data_norm(first_days_per_month,1))

% Change the 'xticklabel' values of the current axes to be your cell array from Problem 3
xticklabels(c_a_3)

% Resize the figure to be wider, so that all the month labels
% on the x-axis are not overlapping.
set(gcf,'units','normalized') % So that it does not depend on the resolution of the user's monitor
set(gcf,'position',[0.1411 0.2843 0.7641 0.3889])
print(fig2,'lab_3_part_5','-depsc','-r600'); % Save figure in a colored .eps format using 600 dpi resolution

%% Part (6)
% Download the file on Canvas called "surprise.txt", which contains
% multiple lines of text, each with a different number of characters.
% First, using the method outlined above, load each line of the file into a cell array in MATLAB.

ct = 0;  % a counter
fid = fopen('surprise.txt', 'r'); % open file for reading
% Keep running a loop until EndOfFile
%   (feof(fid) = 1 if end of the file is reached)
while feof(fid) ~= 1
    ct = ct + 1; % increase counter by 1
    data{ct} = fgetl(fid); % grab entire line of text from file,store as string in the ct^th cell
end
fclose(fid); % close file
%
% ct should be equal to the number of lines of text in the file

% Then, write a code (it shouldn't be long...mine is only 7 lines) that
% will print each line of text into the command window (fprintf), one
% character at a time, and after each full line has been printed, it should
% then print the newline character combination in MATLAB, fprintf('\n'), in
% order to jump to the next line.

for i = 1:ct % Iterate through all lines of text in the file
    for j = 1:length(data{i}) % Iterate through all characters of each line
        fprintf(data{i}(j)) % Print one character at a time
    end
    fprintf('\n') % At this point a full line has been printed, hence the newline character is printed
end

%% Part (7)
% Develop a data set that is 10,000 points long composed of random
% numbers from a normal distribution.  Make a histogram of this data set
% and then overlay a Gaussian curve on top.  You will have to determine how
% to normalize your plot or histogram to get the y-axis to scale appropriately.

rand_data_7 = randn(10000,1); % Create the 10,000 points long vector composed of random numbers from a normal distribution
fig3 = figure(3); % Open a new figure for plotting
pd_7 = fitdist(rand_data_7,'Normal'); % Fit a normal distribution to the generated data
x_values_7 = pd_7.mu-4*pd_7.sigma:0.01:pd_7.mu+4*pd_7.sigma; % Calculate x-coordinates of gaussian curve from -4/+4 standard deviations
y_7 = pdf(pd_7,x_values_7); % Calculate corresponding y-coordinates of fitted curve
yyaxis right % Utilize the right vertical side as a new axis for the PDF
plot(x_values_7,y_7,'Linewidth',2); % Plot the Power Density Function of the fitted normal distribution
set(gca,'YLim',1.5*get(gca,'YLim')) % Increase the limits of the y-axis to avoid overlapping between lines and the legends

yyaxis left % Utilize the left vertical side as a new axis for the histogram
histogram(rand_data_7) % Create a histogram of the generated data, showing the frequency of each element inside the vector
set(gca,'YLim',1.5*get(gca,'YLim')) % Increase the limits of the y-axis to avoid overlapping between lines and the legends

% Create a new legend specifying everything on the figure
legend({' Frequency',[' PDF of Fitted Normal Distribution\newline Mean: ',num2str(pd_7.mu),' - SD: ',num2str(pd_7.sigma)]})
xlabel('Numbers in the Random Dataset'); % Add a label for the x-axis
yyaxis right % Switch back to the right y-axis
ylabel('Normalized Power Density Function'); % Add a label for the right y-axis
yyaxis left % Switch back to the left y-axis
ylabel('Frequency'); % Add a label for the left y-axis
title('LAB \#3 - Part (7) - Histogram - Normal Distribution','FontSize',14,'interpreter','latex') % Create a title for the figure
print(fig3,'lab_3_part_7','-depsc','-r600'); % Save figure in a colored .eps format using 600 dpi resolution

%% Part (8)
% Perform a power spectral density (PSD) calculation using pwelch

t = 0:0.01:10000; % Make a time vector, t, out to 10000 with spacing of 0.01.  
y = 4*cos(2*pi*t)+sin(2*pi*t/0.2)+0.01*randn(1,length(t)); % Generate data
Fs = 100; % Sampling frequency of input vector

% Load Default Values for WINDOW, NOVERLAP, NFFT, based on the information
% from "help pwelch"

% By default, X is divided into the longest possible sections, to get as
% close to but not exceeding 8 segments with 50% overlap.
WINDOW = floor(length(y)/8); % Define a WINDOW size that will lead to 8 segments
NOVERLAP = WINDOW/2; % Set a 50% overlap

% If NFFT is specified as empty, NFFT is set to either 256 or the next
% power of two greater than the length of each section of X, whichever is larger.
NFFT = pow2(floor(log2(WINDOW)+1)); % next power of two greater than the length of each section of X
if NFFT < 256
    % If the next of two greater than the length of each section of y (WINDOW) is smaller than 256, switch to 256
    NFFT = 256;
end

fig4 = figure(4); % Open a new figure for plotting
hold on; % Use "hold on" to plot all three parabolas on the same plot
[Pxx,F] = pwelch(y,WINDOW,NOVERLAP,NFFT,Fs); % Plot the data from pwelch for the default values for the parameters
semilogy(F,Pxx,'Linewidth',2) % Plot on a semilogy plot with Pxx as a function of F (frequency)
xlabel('Frequency (Hz)') % Set a label for the x-axis 
ylabel('Log-PSD [$V^{2}/Hz$]') % Set a label for the y-axis 
title('LAB \#3 - Part (8) - PSD for default parameters of pwelch','FontSize',14,'interpreter','latex') % Create a title for the figure
print(fig4,'lab_3_part_8_default','-depsc','-r600'); % Save figure in a colored .eps format using 600 dpi resolution
%--------------------------------------------------------------------------
% Modify only the paremeter WINDOW
% The smaller the window, the smaller the resolution of the spectrum
[Pxx,F] = pwelch(y,WINDOW,NOVERLAP,NFFT,Fs);
WINDOW_1 = WINDOW/1.5; % Lowest possible value, since NOVERLAP has to be lower than WINDOW
[Pxx_1,F_1] = pwelch(y,WINDOW_1,NOVERLAP,NFFT,Fs);
WINDOW_2 = WINDOW*1.5; % Slightly larger WINDOW than the default value
[Pxx_2,F_2] = pwelch(y,WINDOW_2,NOVERLAP,NFFT,Fs);
WINDOW_3 = WINDOW*5; % Larger WINDOW than the default value
[Pxx_3,F_3] = pwelch(y,WINDOW_3,NOVERLAP,NFFT,Fs);

fig5 = figure(5); % Open a new figure for plotting
hold on; % Use "hold on" to plot all three parabolas on the same plot
semilogy(F,Pxx,'Linewidth',2) % Plot on a semilogy plot with Pxx as a function of F (frequency)
hold on; % Use "hold on" to plot all three parabolas on the same plot
semilogy(F_1,Pxx_1,'Linewidth',2) % Plot on a semilogy plot with Pxx as a function of F (frequency)
hold on; % Use "hold on" to plot all three parabolas on the same plot
semilogy(F_2,Pxx_2,'Linewidth',2) % Plot on a semilogy plot with Pxx as a function of F (frequency)
hold on; % Use "hold on" to plot all three parabolas on the same plot
semilogy(F_3,Pxx_3,'Linewidth',2) % Plot on a semilogy plot with Pxx as a function of F (frequency)
legend('Default','WINDOW/1.5','WINDOW*1.5','WINDOW*5')
title('LAB \#3 - Part (8) - PSD for various values of the WINDOW parameter','FontSize',14,'interpreter','latex') % Create a title for the figure
% Resize the figure to be wider, so that all the month labels
% on the x-axis are not overlapping.
set(gcf,'units','normalized') % So that it does not depend on the resolution of the user's monitor
set(gcf,'position',[0 0.2843 1 0.3889])
xlabel('Frequency (Hz)') % Set a label for the x-axis 
ylabel('Log-PSD [$V^{2}/Hz$]') % Set a label for the y-axis 

% Define second set of axes to zoom over the first peak
ax2 = axes('Position',[0.2 0.6 0.28 0.28],'Box','on'); % Define second set of axes
hold on; % Use "hold on" to plot all three parabolas on the same plot
semilogy(F,Pxx,'Linewidth',2) % Plot on a semilogy plot with Pxx as a function of F (frequency)
hold on; % Use "hold on" to plot all three parabolas on the same plot
semilogy(F_1,Pxx_1,'Linewidth',2) % Plot on a semilogy plot with Pxx as a function of F (frequency)
hold on; % Use "hold on" to plot all three parabolas on the same plot
semilogy(F_2,Pxx_2,'Linewidth',2) % Plot on a semilogy plot with Pxx as a function of F (frequency)
hold on; % Use "hold on" to plot all three parabolas on the same plot
semilogy(F_3,Pxx_3,'Linewidth',2) % Plot on a semilogy plot with Pxx as a function of F (frequency)
axis(ax2,[0.995 1.005 ax2.YLim]); % Define the axis limits of the second set of axes to observe a close neighborhood of the target

% Define third set of axes to zoom over the second peak
ax3 = axes('Position',[0.3 0.25 0.28 0.28],'Box','on'); % Define second set of axes
hold on; % Use "hold on" to plot all three parabolas on the same plot
semilogy(F,Pxx,'Linewidth',2) % Plot on a semilogy plot with Pxx as a function of F (frequency)
hold on; % Use "hold on" to plot all three parabolas on the same plot
semilogy(F_1,Pxx_1,'Linewidth',2) % Plot on a semilogy plot with Pxx as a function of F (frequency)
hold on; % Use "hold on" to plot all three parabolas on the same plot
semilogy(F_2,Pxx_2,'Linewidth',2) % Plot on a semilogy plot with Pxx as a function of F (frequency)
hold on; % Use "hold on" to plot all three parabolas on the same plot
semilogy(F_3,Pxx_3,'Linewidth',2) % Plot on a semilogy plot with Pxx as a function of F (frequency)
axis(ax3,[4.995 5.005 ax3.YLim(1) 500]); % Define the axis limits of the second set of axes to observe a close neighborhood of the target
print(fig5,'lab_3_part_8_window','-depsc','-r600'); % Save figure in a colored .eps format using 600 dpi resolution
%--------------------------------------------------------------------------
% Modify only the paremeter NOVERLAP
[Pxx,F] = pwelch(y,WINDOW,NOVERLAP,NFFT,Fs);
NOVERLAP_1 = NOVERLAP/2;
[Pxx_1,F_1] = pwelch(y,WINDOW,NOVERLAP_1,NFFT,Fs);
NOVERLAP_2 = NOVERLAP/4;
[Pxx_2,F_2] = pwelch(y,WINDOW,NOVERLAP_2,NFFT,Fs);
NOVERLAP_3 = NOVERLAP/8;
[Pxx_3,F_3] = pwelch(y,WINDOW,NOVERLAP_3,NFFT,Fs);
NOVERLAP_4 = NOVERLAP*1.5 % Highest possible value is WINDOW-1, since NOVERLAP has to be lower than WINDOW
[Pxx_4,F_4] = pwelch(y,WINDOW,NOVERLAP_4,NFFT,Fs);

fig6 = figure(6); % Modify only the 'WINDOW' parameter
hold on; % Use "hold on" to plot all three parabolas on the same plot
semilogy(F,Pxx,'Linewidth',2)
hold on; % Use "hold on" to plot all three parabolas on the same plot
semilogy(F_1,Pxx_1,'Linewidth',2)
hold on; % Use "hold on" to plot all three parabolas on the same plot
semilogy(F_2,Pxx_2,'Linewidth',2)
hold on; % Use "hold on" to plot all three parabolas on the same plot
semilogy(F_3,Pxx_3,'Linewidth',2)
hold on; % Use "hold on" to plot all three parabolas on the same plot
semilogy(F_4,Pxx_4,'Linewidth',2)
legend('Default','NOVERLAP/2','NOVERLAP/4','NOVERLAP/8','NOVERLAP*1.5')
title('LAB \#3 - Part (8) - PSD for various values of the NOVERLAP parameter','FontSize',14,'interpreter','latex') % Create a title for the figure
% Resize the figure to be wider, so that all the month labels
% on the x-axis are not overlapping.
set(gcf,'units','normalized') % So that it does not depend on the resolution of the user's monitor
set(gcf,'position',[0 0.2843 1 0.3889])
xlabel('Frequency (Hz)') % Set a label for the x-axis 
ylabel('Log-PSD [$V^{2}/Hz$]') % Set a label for the y-axis 

% Define second set of axes to zoom over the first peak
ax2 = axes('Position',[0.2 0.6 0.28 0.28],'Box','on'); % Define second set of axes
hold on; % Use "hold on" to plot all three parabolas on the same plot
semilogy(F,Pxx,'Linewidth',2)
hold on; % Use "hold on" to plot all three parabolas on the same plot
semilogy(F_1,Pxx_1,'Linewidth',2)
hold on; % Use "hold on" to plot all three parabolas on the same plot
semilogy(F_2,Pxx_2,'Linewidth',2)
hold on; % Use "hold on" to plot all three parabolas on the same plot
semilogy(F_3,Pxx_3,'Linewidth',2)
hold on; % Use "hold on" to plot all three parabolas on the same plot
semilogy(F_4,Pxx_4,'Linewidth',2)
axis(ax2,[0.995 1.005 ax2.YLim]); % Define the axis limits of the second set of axes to observe a close neighborhood of the target

% Define third set of axes to zoom over the second peak
ax3 = axes('Position',[0.3 0.25 0.28 0.28],'Box','on'); % Define second set of axes
hold on; % Use "hold on" to plot all three parabolas on the same plot
semilogy(F,Pxx,'Linewidth',2)
hold on; % Use "hold on" to plot all three parabolas on the same plot
semilogy(F_1,Pxx_1,'Linewidth',2)
hold on; % Use "hold on" to plot all three parabolas on the same plot
semilogy(F_2,Pxx_2,'Linewidth',2)
hold on; % Use "hold on" to plot all three parabolas on the same plot
semilogy(F_3,Pxx_3,'Linewidth',2)
hold on; % Use "hold on" to plot all three parabolas on the same plot
semilogy(F_4,Pxx_4,'Linewidth',2)
axis(ax3,[4.995 5.005 ax3.YLim(1) 500]); % Define the axis limits of the second set of axes to observe a close neighborhood of the target
print(fig6,'lab_3_part_8_noverlap','-depsc','-r600'); % Save figure in a colored .eps format using 600 dpi resolution
%--------------------------------------------------------------------------
% Modify only the paremeter NFFT
[Pxx,F] = pwelch(y,WINDOW,NOVERLAP,NFFT,Fs);
NFFT_1 = round(NFFT/2.5);
[Pxx_1,F_1] = pwelch(y,WINDOW,NOVERLAP,NFFT_1,Fs);
% NFFT_2 = round(NFFT/10);
% [Pxx_2,F_2] = pwelch(y,WINDOW,NOVERLAP,NFFT_2,Fs);
NFFT_3 = NFFT*10 % Highest possible value is WINDOW-1, since NOVERLAP has to be lower than WINDOW
[Pxx_3,F_3] = pwelch(y,WINDOW,NOVERLAP,NFFT_3,Fs);

fig7 = figure(7) % Modify only the 'WINDOW' parameter
hold on; % Use "hold on" to plot all three parabolas on the same plot
semilogy(F,Pxx,'Linewidth',2)
hold on; % Use "hold on" to plot all three parabolas on the same plot
semilogy(F_1,Pxx_1,'Linewidth',2)
% hold on; % Use "hold on" to plot all three parabolas on the same plot
% semilogy(F_2,Pxx_2,'Linewidth',2)
hold on; % Use "hold on" to plot all three parabolas on the same plot
semilogy(F_3,Pxx_3,'Linewidth',2)
legend('Default','NFF5/2.5','NFT*5','NFT*5')
title('LAB \#3 - Part (8) - PSD for various values of the NFFT parameter','FontSize',14,'interpreter','latex') % Create a title for the figure
% Resize the figure to be wider, so that all the month labels
% on the x-axis are not overlapping.
set(gcf,'units','normalized') % So that it does not depend on the resolution of the user's monitor
set(gcf,'position',[0 0.2843 1 0.3889])
xlabel('Frequency (Hz)') % Set a label for the x-axis 
ylabel('Log-PSD [$V^{2}/Hz$]') % Set a label for the y-axis 

% Define second set of axes to zoom over the first peak
ax2 = axes('Position',[0.2 0.6 0.28 0.28],'Box','on'); % Define second set of axes
hold on; % Use "hold on" to plot all three parabolas on the same plot
semilogy(F,Pxx,'Linewidth',2) % Plot on a semilogy plot with Pxx as a function of F (frequency)
hold on; % Use "hold on" to plot all three parabolas on the same plot
semilogy(F_1,Pxx_1,'Linewidth',2) % Plot on a semilogy plot with Pxx as a function of F (frequency)
% hold on; % Use "hold on" to plot all three parabolas on the same plot
% semilogy(F_2,Pxx_2,'Linewidth',2) % Plot on a semilogy plot with Pxx as a function of F (frequency)
hold on; % Use "hold on" to plot all three parabolas on the same plot
semilogy(F_3,Pxx_3,'Linewidth',2) % Plot on a semilogy plot with Pxx as a function of F (frequency)
axis(ax2,[0.995 1.005 ax2.YLim]); % Define the axis limits of the second set of axes to observe a close neighborhood of the target

% Define third set of axes to zoom over the second peak
ax3 = axes('Position',[0.3 0.25 0.28 0.28],'Box','on'); % Define second set of axes
hold on; % Use "hold on" to plot all three parabolas on the same plot
semilogy(F,Pxx,'Linewidth',2) % Plot on a semilogy plot with Pxx as a function of F (frequency)
hold on; % Use "hold on" to plot all three parabolas on the same plot
semilogy(F_1,Pxx_1,'Linewidth',2) % Plot on a semilogy plot with Pxx as a function of F (frequency)
% hold on; % Use "hold on" to plot all three parabolas on the same plot
% semilogy(F_2,Pxx_2,'Linewidth',2) % Plot on a semilogy plot with Pxx as a function of F (frequency)
hold on; % Use "hold on" to plot all three parabolas on the same plot
semilogy(F_3,Pxx_3,'Linewidth',2) % Plot on a semilogy plot with Pxx as a function of F (frequency)
axis(ax3,[4.995 5.005 ax3.YLim(1) 500]); % Define the axis limits of the second set of axes to observe a close neighborhood of the target
print(fig7,'lab_3_part_8_nfft','-depsc','-r600'); % Save figure in a colored .eps format using 600 dpi resolution
