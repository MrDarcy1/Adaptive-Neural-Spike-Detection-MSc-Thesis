load('realDataWithLFP_1.mat')
load('spike_location_1.mat')
N=50000;
data=round((1e7/(4))*data(1:N));
window_length = 16;
Mean=filter(ones(1,window_length),1,data);
Mean=Mean(window_length/2+1:end)/window_length;


mean_buffer_mean = round(mean(data(1:window_length)));
MEAN = mean_buffer_mean;

mean_buffer = data(1:window_length);
mean_buffer_end = 1;
% countMax = 500;
count = 0;
var = [];
for countMax = 100:5:200
    MEAN = MEAN(1);
    for i = 17:N
        temp = mean_buffer(mean_buffer_end);
        mean_buffer(mean_buffer_end) = data(i);
        if (count == countMax)
            count = 0;
            mean_buffer_mean = round(mean(mean_buffer));
        else
            mean_buffer_mean = mean_buffer_mean - round((temp - mean_buffer(mean_buffer_end))/16);
            count = count + 1;
        end
        if mean_buffer_end == window_length
            mean_buffer_end = 1;
        else
            mean_buffer_end = mean_buffer_end + 1;
        end
        MEAN = [MEAN mean_buffer_mean];
    end
    var = [var sum((Mean(8:end)'-MEAN).^2)/length(MEAN)];
end
% 
% figure(1);
% plot(Mean(8:end));hold on;
% plot(MEAN);hold off
    