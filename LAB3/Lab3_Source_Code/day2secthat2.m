function [num_sec] = day2secthat2(mat_dates_v)
% function [num_sec] = day2secthat2(mat_dates_v)

% Description: This function receives as an input a vector of MATLAB dates
% and returns as an output a vector containing the number of seconds that
% have elapsed since the 1st value (date) in the input vector. This a
% second simpler version of the "day2secthat" function

% Calculate elapsed days of all elements wrt to the first date in the input
days_elapsed = mat_dates_v-mat_dates_v(1);

% Calculate elapsed hours of all elements wrt to the first date in the input
hours_elapsed = days_elapsed*24; % Each day has 24 hours

% Calculate elapsed minutes of all elements wrt to the first date in the input
minutes_elapsed = hours_elapsed*60; % Each hour has 60 minutes

% Calculate elapsed seconds of all elements wrt to the first date in the input
num_sec = minutes_elapsed*60; % Each minute has 60 seconds
end