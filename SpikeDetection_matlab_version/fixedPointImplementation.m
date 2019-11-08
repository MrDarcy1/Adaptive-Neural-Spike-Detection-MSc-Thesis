% Important parameters: c, c_init, updateFreq
% 100000 samples are taken. ~4s. 
%%
addpath('./data')
addpath('./functions')
addpath('./utils')

%% Load Data
% data = signal;
% spike_location = corrSpikes;
load('realDataWithLFP_3.mat')
load('spike_location_3.mat')
N=100000;
data=round((1e7/(4))*data(1:N));
spike_location=spike_location(spike_location<=N);


% figure(2)
% plotSpikes(spike_location,data)
% legend('Introcellalur Recording','WaveClus Detected Spikes')
% xlabel('Time/s')
% ylabel('Amplitude')
% spike_location=spike_location(spike_location>N)-N;
fs=24414;
%% Add Noise
load('noise_base.mat');
SNR=5;
lambda=20; %spike freq
cells=3;
[noise_data,noise]=addNoise(data,SNR,12);

% [noise_data,noise,backgroundActNum,backgroundActLoc] = addNoisePossion(data,noise_base,SNR,lambda,cells,fs);
%     data = noise_data/100;
% %     data=1e4*double(data(1:N));
%     fid=fopen(['noise_data.txt'],'wt'); %写的方式打开文件（若不存在，建立文件）；
%     fprintf(fid,'%f\n',data);  % %d 表示以整数形式写入数据，这正是我想要的；
% figure(1)
% subplot(3,1,1)
% plotSpikes(spike_location,data)
% title('Original Data')
% xlabel('Time Steps')
% ylabel('Amplitude')
% % ylim([-3e-4,6e-4])
% xlim([1,5000])
% subplot(2,1,1)
% plotSpikes(spike_location,noise_data)
% ylim([-0.4e-3,0.6e-3]);
% xlim([1,5000])
% % ylim([-3e-4,6e-4])
% title('Signal with 5dB Gaussian Noise')
% xlabel('Time/s')
% ylabel('Amplitude')
% subplot(2,1,2)
% plot(0:1/fs:(N-1)/fs,noise)
% xlim([0,N/fs])
% ylim([-3e-4,6e-4])
% title(['White Gaussian Noise - SNR: ',num2str(SNR),'dB'])
% xlabel('Time/s')
% ylabel('Amplitude')
% ylim([-0.5e-3,0.5e-3]);

%% Main


% data = round(noise_data);


countFreq = 1000;
updateFreq = 15000;
c= 24;
c_init = 8;
%out
ASO = [];
THR = [];
DEMEAN = [];
MEAN = [];
DATA_FOR_THR = [];
spikes_detected = [];
interval = [];
%parameters
count = 0;
adding = 0;
Sum = 0;
mean_buffer_length = 16;
mean_buffer_end = 1;
previous_demean = 0;
middle_value = 0;
preprevious_demean = 0;
ppreprevious_demean = 0;
pppreprevious_demean = 0;
mean_buffer = zeros(mean_buffer_length,1);

aso = 0;

thr_buffer_length = 64;
thr_buffer_end = 1;
thr_buffer = zeros(thr_buffer_length,1);

update_time = updateFreq;
detected_time = 30;
hold_time = 20;


detected = detected_time;
update = update_time;
hold_data = hold_time;
% 
after = 15;

%fill buffers
SUM = 0;
for i = 1:mean_buffer_length
    mean_buffer(mean_buffer_end) = data(i);
    SUM = SUM + mean_buffer(mean_buffer_end);
    mean_buffer_end = mean_buffer_end + 1;
end

mean_buffer_end = 1;
mean_buffer_mean = round(SUM/mean_buffer_length);

SUM = 0;
for i = 1:thr_buffer_length
    temp = mean_buffer(mean_buffer_end);
    mean_buffer(mean_buffer_end) = data(i+16);
    demean = mean_buffer(mean_buffer_end) - mean_buffer_mean;
    mean_buffer_mean = mean_buffer_mean - round((temp - mean_buffer(mean_buffer_end))/16);
    if mean_buffer_end == mean_buffer_length
        mean_buffer_end = 1;
    else
        mean_buffer_end = mean_buffer_end + 1;
    end
    aso = round(pow2(nextpow2(demean)-1)*abs(demean-preprevious_demean)/2^7);
    ppreprevious_demean = preprevious_demean;
    preprevious_demean=previous_demean;
    previous_demean = demean;

    thr_buffer(thr_buffer_end) = aso;
    thr_buffer_end = thr_buffer_end + 1;
end
threshold = c_init*median(thr_buffer);
sampleNum = N;


% processing
for i = 1:sampleNum
    temp = mean_buffer(mean_buffer_end);
    mean_buffer(mean_buffer_end) = data(i);
    if (count == countFreq)
        count = 0;
        mean_buffer_mean = round(mean(mean_buffer));
    else
        mean_buffer_mean = (mean_buffer_mean - round((temp - mean_buffer(mean_buffer_end))/16));
        count = count + 1;
    end
    demean = mean_buffer(mean_buffer_end) - mean_buffer_mean;
    if mean_buffer_end == mean_buffer_length
        mean_buffer_end = 1;
    else
        mean_buffer_end = mean_buffer_end + 1;
    end
    aso = round(pow2(nextpow2(demean)-1)*abs(demean-preprevious_demean)/2^7);
    
    ppreprevious_demean = preprevious_demean;
    preprevious_demean=previous_demean;
    previous_demean = demean;
    
    if(detected == detected_time)
        if(aso>threshold) %detection

            spikes_detected=[spikes_detected i]; %record spike location
            interval=[interval;[i,i+after]]; % record spike interval
            detected = 0;
        end
        
        if (update < update_time)
            update = update + 1;
        elseif (adding < thr_buffer_length)
            if(aso<=threshold / 2)
                Sum = Sum + aso;
                adding = adding + 1;
            end
        else
            if(aso <= threshold / 2)
                Sum = Sum + aso;
                threshold = round(Sum*c/thr_buffer_length);
                Sum = 0;
                adding = 0;
                update = 0;
            end
        end          
    else
        detected = detected + 1;
        update = update + 1;
    end
    MEAN = [MEAN, mean_buffer_mean];
    DEMEAN = [ DEMEAN, demean];
    ASO = [ASO, aso];
    THR = [THR, threshold];
end

% figure(1);
% plot(ASO);hold on
% plot(DATA_FOR_THR(64:end));hold off
% 
% figure(2)
% plot(ASO(65:end));hold on
% plot(THR);hold off
% 
% figure(3);
% plot(data);hold on
% plot(MEAN(8:end));hold off;
% 
% figure(4)
% plotSpikes(spikes_detected,ASO);hold on
% plot(THR);hold off


[FP,FN,TP]=locationCompare(spike_location,interval,spikes_detected);

% Sens(l,j,k)= length(TP)/(length(TP)+length(FN)); % found is correct
% FDR(l,j,k) = length(FP)/(length(FP)+length(TP));% not find
% Acc(l,j,k) = length(TP)/(length(TP)+length(FN)+length(FP));
Sens1= length(TP)/(length(TP)+length(FN)) % found is correct
FDR1 = length(FP)/(length(FP)+length(TP))% not find
Acc1 = length(TP)/(length(TP)+length(FN)+length(FP))
% Sens(l,j,k)
% FDR(l,j,k)
% Acc(l,j,k)
%     k = k + 1
% 
% figure(5)
visualisation(ASO,THR,FP,FN,TP,0)

%         end
%         j = j + 1
%     end
%     l = l + 1
% end
% figure(5)
% 
  
