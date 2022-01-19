function [num_sec] = day2secthat(mat_dates_v)
% function [num_sec] = day2secthat(mat_dates_v)

% Description: This function receives as an input a vector of MATLAB dates
% and returns as an output a vector containing the number of seconds that
% have elapsed since the 1st value (date) in the input vector

mat2str = datestr(mat_dates_v); % Convert the input vector to string format
str2vec = datevec(mat2str); % Convert the string format to a vector of components
num_sec = []; % Initialize the output vector that will contain the elapsed number of seconds
for i = 1:size(str2vec,1) % Iterate through all dates (elements) of the input vector
    % Calculate the elapsed number of seconds of each date wrt to the
    % first one and add it to the output vector.
    num_sec = [num_sec etime(str2vec(i,:),str2vec(1,:))];
end
end