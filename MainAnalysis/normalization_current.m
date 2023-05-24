function [data_norm] = normalization_current(data, spotID, max_val)
% normalizes data between 0 and 1 by subtracting the baseline and dividing
% by the difference between the maximum and baseline
% will only normalize if the difference between the baseline and maximum is
% above a specified threshold

baseline = mean(data(3:7));

% if fluorescence diff between start and max is large enough, 
% normalize to 0-1
% otherwise just normalize by baseline
if strcmp(spotID, 'NC')
    data_norm = data ./ baseline - 1;
elseif strcmp(spotID, 'PC')
    data_norm = data ./ baseline;
else
    data_norm = (data - baseline) / max_val;
end
    
end

