function [demean_data,Mean] = extractMean(data,window_length)
    Mean=filter(ones(1,window_length),1,data);
    Mean=Mean(window_length/2+1:end);
    demean_data=data((1:end-window_length/2))-Mean/window_length;
end

