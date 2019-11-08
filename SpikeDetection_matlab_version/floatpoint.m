% Important parameters: c, c_init, updateFreq
% 100000 samples are taken. ~4s. 
%%

addpath('./data')
addpath('./functions')
addpath('./utils')

%% Load Data
load('realDataWithLFP_3.mat')
load('spike_location_3.mat')
N=100000;
data=1e7*data(1:N);
% [b,a] = butter(2,300/3500,'high');
% data = filter(b,a,data);
% data = data(1:N);
% subplot(4,1,1)
spike_location=spike_location(spike_location<=N);

% plotSpikes(data,spike_location)
% xlabel('Time/s')
% title('Introcellular Signal')
% spike_location=spike_location(spike_location>N)-N;
fs=24414;
%% Add Noise
% load('noise_base.mat');
% SNR=inf;
% lambda=20; %spike freq
% cells=3;
[noise_data,noise]=addNoise(data,3,-5);
% 
% % [noise_data,noise,backgroundActNum,backgroundActLoc] = addNoisePossion(data,noise_base,SNR,lambda,cells,fs);
% %     data = noise_data/100;
% % %     data=1e4*double(data(1:N));
% %     fid=fopen(['noise_data.txt'],'wt');
% %     fprintf(fid,'%f\n',data);  
% figure(1)
% subplot(3,1,1)
% plotSpikes(spike_location,data)
% title('Original Data')
% xlabel('Time Steps')
% ylabel('Amplitude')
% % ylim([-3e-4,6e-4])
% xlim([0,(N-1)/fs])
% subplot(3,1,2)
% plotSpikes(spike_location,noise_data)
% xlim([0,(N-1)/fs])
% % ylim([-3e-4,6e-4])
% title('Signal with Noise')
% xlabel('Time Steps')
% ylabel('Amplitude')
% subplot(3,1,3)
% plot(0:1/fs:(N-1)/fs,noise)
% xlim([0,(N-1)/fs])
% % ylim([-3e-4,6e-4])
% title(['Background Activity - SNR: ',num2str(SNR),'dB'])
% xlabel('Time Steps')
% ylabel('Amplitude')
%% Main
% data = noise_data;


countFreq = 1000;
updateFreq = 15000;
c = 24;
%out
ASO = [];
THR1 = [];
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

thr_buffer_length = 65;
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
mean_buffer_mean = (SUM/mean_buffer_length);

SUM = 0;
for i = 1:thr_buffer_length
    temp = mean_buffer(mean_buffer_end);
    mean_buffer(mean_buffer_end) = data(i+16);
    demean = mean_buffer(mean_buffer_end) - mean_buffer_mean;
    mean_buffer_mean = mean_buffer_mean - ((temp - mean_buffer(mean_buffer_end))/16);
    if mean_buffer_end == mean_buffer_length
        mean_buffer_end = 1;
    else
        mean_buffer_end = mean_buffer_end + 1;
    end
%     demean = data(i);
    aso = abs(demean*(demean-preprevious_demean)/2^7);
%     aso = abs(demean^2 - previous_demean * (data(i+17)-mean_buffer_mean)/2^7) ;

    ppreprevious_demean = preprevious_demean;
    preprevious_demean=previous_demean;
    previous_demean = demean;

    thr_buffer(thr_buffer_end) = aso;
    thr_buffer_end = thr_buffer_end + 1;
end
threshold = 15*median(thr_buffer);
sampleNum = N-2;
% processing
for i = 1:sampleNum
    temp = mean_buffer(mean_buffer_end);
    mean_buffer(mean_buffer_end) = data(i);
    if (count == countFreq)
        count = 0;
        mean_buffer_mean = (mean(mean_buffer));
    else
        mean_buffer_mean = (mean_buffer_mean - ((temp - mean_buffer(mean_buffer_end))/16));
        count = count + 1;
    end
    demean = abs(mean_buffer(mean_buffer_end) - mean_buffer_mean);
% demean = data(i);
    if mean_buffer_end == mean_buffer_length
        mean_buffer_end = 1;
    else
        mean_buffer_end = mean_buffer_end + 1;
    end
    aso = abs(demean*(demean-preprevious_demean)/2^7);
    
%     aso = abs(demean^2 - preprevious_demean * (data(i+2)-mean_buffer_mean)/2^7) ;
    ppreprevious_demean = preprevious_demean;
    preprevious_demean=previous_demean;
    previous_demean = demean;
    
    if(detected == detected_time)
        if(aso>threshold) %detection
%             if i+after<sampleNum
%                 [~,idx]=max(aso(i:i+after));
%             else %last few data points
%                 [~,idx]=max(aso(i:sampleNum));
%             end
            spikes_detected=[spikes_detected i]; %record spike location
            interval=[interval;[i,i+after]]; % record spike interval
            detected = 0;
        end
        
        if (update < update_time)
            DATA_FOR_THR = [DATA_FOR_THR, 0];

            update = update + 1;
        elseif (adding < thr_buffer_length)
            if(aso<=threshold / 2)
                DATA_FOR_THR = [DATA_FOR_THR, aso];
                Sum = Sum + aso;
                adding = adding + 1;
            end
        else
             

            if(aso <= threshold / 2)
                DATA_FOR_THR = [DATA_FOR_THR, aso];
                Sum = Sum + aso;
                threshold = (Sum*c/thr_buffer_length);
                Sum = 0;
                adding = 0;
                update = 0;
            end
        end          
    else
        DATA_FOR_THR = [DATA_FOR_THR, 0];

        detected = detected + 1;
        update = update + 1;
    end
    MEAN = [MEAN, mean_buffer_mean];
    DEMEAN = [ DEMEAN, demean];
    ASO = [ASO, aso];
    THR1 = [THR1, threshold];
end
% 

% figure(1);
% plot(0:1/fs:(length(ASO)-1)/fs,ASO,'Color','g','LineWidth',0.5);hold on
% % plot(0:1/fs:(length(DATA_FOR_THR1)-1)/fs,DATA_FOR_THR1,'Color','b','LineWidth',0.5);hold on
% plot(0:1/fs:(length(DATA_FOR_THR)-1)/fs,DATA_FOR_THR,'Color','r','LineWidth',0.5);hold off
% xlim([0,(length(DATA_FOR_THR)-1)/fs])
% xlabel('Time/s')
% % legend('ASO Processed Data','Spike Eliminted','Spike Eliminated and Sub-thresholded')
% title('Data used for thresholding')


% 
% figure(10)
% plotSpikes(spikes_detected,ASO);hold on
% plot(0:1/fs:(length(ASO)-1)/fs,ASO);hold on
% plot(0:1/fs:(length(THR1)-1)/fs,THR1);hold off
% plot(0:1/fs:(length(THR)-1)/fs,THR);hold off
% xlim([0,(length(THR)-1)/fs])
% xlabel('Time/s')
% legend('ASO Processed Data','Threshold ','Spike Eliminated and Sub-thresholded Threshold')
% title('Thresholds')

% figure(3);
% plot(data);hold on
% plot(MEAN(8:end));hold off;
% 
% figure(4)
% subplot(4,1,1)
% plotSpikes(spike_location,data);
% legend('Intracellular Signal','Spike Locations')
% xlabel('Time/s')
% xlim([0,(length(THR1)-1)/fs])
% title('Intracellular Signal - SNR decreases from 3dB to -5dB')
% 
% subplot(4,1,2)
% plotSpikes(spike_location,DEMEAN);
% legend('Local Field Potential Removed Signal','Spike Locations')
% xlabel('Time/s')
% xlim([0,(length(THR1)-1)/fs])
% title('Local Field Potential Removed Signal')
% 
% subplot(4,1,3)
% plotSpikes(spike_location,ASO);
% legend('Spike-Emphasised Signal','Spike Locations')
% xlabel('Time/s')
% xlim([0,(length(THR1)-1)/fs])
% title('Spike Emphasised Signal')
% 
% subplot(4,1,4)
% plotSpikes(spikes_detected,ASO);hold on
% plot(0:1/fs:(length(THR1)-1)/fs,THR1,'Color','r','LineWidth',1);hold off
% xlim([0,(length(THR1)-1)/fs])
% xlabel('Time/s')
% legend('Emphasised Signal','Detected Spikes','Threshold')%,'Spike Eliminated and Sub-thresholded Threshold')
% title('Statistic Threshold')

[FP,FN,TP]=locationCompare(spike_location,interval,spikes_detected);

% Sens(l)= length(TP)/(length(TP)+length(FN)) % found is correct
% FDR(l) = length(FP)/(length(FP)+length(TP))% not find
% Acc(l) = length(TP)/(length(TP)+length(FN)+length(FP))
Sens= length(TP)/(length(TP)+length(FN)) % found is correct
FDR = length(FP)/(length(FP)+length(TP))% not find
Acc = length(TP)/(length(TP)+length(FN)+length(FP))
% Sens(l,j,k)
% FDR(l,j,k)
% Acc(l,j,k)
%     k = k + 1
% 
% figure(5)
visualisation(ASO,THR1,FP,FN,TP,0)
% xlim([0,(length(THR1)-1)/fs])

%         end
%         j = j + 1
%     end
%     l = l + 1
% end
% figure(5)
% 
%%
% for i = 1:67
% plot(interval(i,1):interval(i,2),ASO(interval(i,1):interval(i,2)),'Color',COLORS(1,:))
% hold on
% end
% hold off

