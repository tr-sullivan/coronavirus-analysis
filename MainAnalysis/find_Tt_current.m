function [Tt, threshold_val] = find_Tt_current(Data_ct, time_ct, y_max, threshold, spotID)
% finds the location at which the time series crosses a user-defined
% threshold using linear interpolation

if strcmp(spotID, 'NC') || strcmp(spotID, 'PC')
    Tt = NaN;
    threshold_val = NaN;
    return
end

% set threshold as fraction of range from baseline
min_val= mean(Data_ct(2:5));
range = y_max - min_val;
threshold_val = min_val + range*threshold/100;

II = find(Data_ct>threshold_val);
I = find(Data_ct<threshold_val);

if isempty(II) || length(II) == 1 || isempty(I) | II(end) < I(1)
    Tt = NaN;
    threshold_val = NaN;
else
    for i = length(II):-1:2
        if II(i-1) ~= II(i)-1
            break;
        end
        i = 1;
    end
    val_2 = II(i);
    val_1 = val_2 - 1;
    data_x_int(1,1) = Data_ct (val_1);
    data_x_int(2,1) = Data_ct (val_2);
    data_y_int(1,1) = time_ct (val_1);
    data_y_int(2,1) = time_ct (val_2);
    Tt= interp1(data_x_int, data_y_int, threshold_val)/60;
    Tt = round(Tt,2);
end

end

