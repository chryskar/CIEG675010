function [min_ind,max_ind] = min_max(y)
% function [min_ind,max_ind] = min_max(y)

% Description: This function receives a smooth vector of data as an input.
% Then the first derivative test is utilized to calculate the local extrema
% of the vector. Finally, the function returns the indices of the input
% vector where these maxima and minima occur.

% Since the input vector of data is smooth, we can differentiate it
y_div = diff(y); %Since we not dividing by time (dt), it is not a derivative, but this vector has the same sign as the derivative

min_ind = []; % Initialize a vector that will contain the indices of the local minima
max_ind = []; % Initialize a vector that will contain the indices of the local maxima
for i = 2:length(y_div) % Iterate through all elements of the "derivative" vector
    if y_div(i) > 0 && y_div(i-1) < 0       % A sign change from positive to negative in the "derivative" indicates a local minimum
        min_ind = [min_ind i];              % Append the vector with all local minima with the new one
    elseif y_div(i) < 0 && y_div(i-1) > 0   % A sign change from negative to positive in the "derivative" indicates a local maximum
        max_ind = [max_ind i];              % Append the vector with all local maxima with the new one
    end
end
end