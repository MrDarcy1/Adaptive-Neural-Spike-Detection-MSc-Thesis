% parameters to test. updatefreq, countFreq, parameters, sampling level,
% sampling freq

addpath('./data')
addpath('./functions')
addpath('./utils')

%% Load Data
for i = 1:3

% plotSpikes(data,spike_location)
% spike_location=spike_location(spike_location>N)-N;
fs=24414;
end
%% Add Noise
load('noise_base.mat');
lambda=20; %spike freq
cells=3;
N=100000;
% [noise_data,noise]=addNoise(data,noise_base,SNR);

j = 1;
for dataNum = 1:3
    load(['realDataWithLFP_',num2str(j),'.mat'])
    load(['spike_location_',num2str(j),'.mat'])
    Data=round((1e7/(4))*data(1:N));
    % data = data(1:N);
    spike_location=spike_location(spike_location<=N);
%     spike_location=spike_location(spike_location>3000)-3000;
    k = 1;
    for SNR = 5:0.5:20.5
        [data,noise,backgroundActNum,backgroundActLoc] = addNoisePossion(Data,noise_base,SNR,lambda,cells,fs);


countFreq = 1000;
updateFreq = 800;
c= 27;
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
detected_time = 20;
hold_time = 20;


detected = detected_time;
update = update_time;
hold_data = hold_time;
% 
after = 10;

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
    SUM = SUM + aso;
    thr_buffer_end = thr_buffer_end + 1;
end
thr_buffer_end = 1;
thr_buffer_mean = round(SUM / thr_buffer_length);
threshold = thr_buffer_mean;
sampleNum = N;


% processing
for i = 1:sampleNum
    temp = mean_buffer(mean_buffer_end);
    mean_buffer(mean_buffer_end) = data(i);
    if (count == countFreq)
        count = 0;
        mean_buffer_mean = round(mean(mean_buffer));
        thr_buffer_mean = round(mean(thr_buffer));
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
        if(aso <= threshold) %not detected
            if(aso > round(threshold / 2))
                hold_data = 0;
            end
            if(update == update_time && hold_data == hold_time )
                update = 0;
                threshold = c*thr_buffer_mean;
            else
                if (update ~= update_time)
                    update = update + 1;
                end
                if(hold_data ~= hold_time)
                    hold_data = hold_data+1;     
                end
            end
        else %detected
            if i+after<sampleNum
                [~,idx]=max(data(i:i+after));
            else %last few data points
                [~,idx]=max(data(i:sampleNum));
            end
            spikes_detected=[spikes_detected i+idx-1]; %record spike location
            interval=[interval;[i,i+after]]; % record spike interval
            detected = 0;
            hold_data = 0;
            update = update_time;
        end
    else
        detected = detected + 1;
        hold_data = hold_data + 1;
    end
    if(hold_data == hold_time)
        thr_buffer_mean = abs(round((thr_buffer_mean*64-thr_buffer(thr_buffer_end)+abs(aso))/64));
        thr_buffer(thr_buffer_end) = aso;
    else
        temp = thr_buffer_mean;
        thr_buffer_mean = abs(round((thr_buffer_mean*64-thr_buffer(thr_buffer_end)+thr_buffer_mean)/64));
        thr_buffer(thr_buffer_end) = temp;
    end
    if(thr_buffer_mean == 0)
        thr_buffer_mean = 2;
    end
    if (thr_buffer_end == thr_buffer_length)
        thr_buffer_end = 1;
    else
        thr_buffer_end = thr_buffer_end + 1;
    end
    MEAN = [MEAN, mean_buffer_mean];
    DATA_FOR_THR = [DATA_FOR_THR, thr_buffer(thr_buffer_end)];
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
% plotSpikes(spikes_detected,ASO)


[FP,FN,TP]=locationCompare(spike_location,interval,spikes_detected);

% Sens(l,j,k)= length(TP)/(length(TP)+length(FN)); % found is correct
% FDR(l,j,k) = length(FP)/(length(FP)+length(TP));% not find
% Acc(l,j,k) = length(TP)/(length(TP)+length(FN)+length(FP));
Sens(j,k)= length(TP)/(length(TP)+length(FN)) % found is correct
FDR(j,k) = length(FP)/(length(FP)+length(TP))% not find
Acc(j,k) = length(TP)/(length(TP)+length(FN)+length(FP))
Sens(j,k)
FDR(j,k)
Acc(j,k)
figure(1)
visualisation(ASO,THR,FP,FN,TP,0)
k =k + 1;
end
j = j + 1;
end

%     k = k + 1
% 

%         end
%         j = j + 1
%     end
%     l = l + 1
% end
% figure(5)
% 
%%
figure(1)
for i = 1:3
    subplot(1,3,i)
    plot(5:0.5:20.5,smooth(Acc(i,:)));hold on;grid on;
    plot(5:0.5:20.5,smooth(fpAcc(i,:)));hold off
    title(['data ', num2str(i),' Accuracy']);
    legend('Fixed Point Result','Float Point Result')
    xlabel('SNR')
    ylim([0,1])
end
saveas(1,'./plt/fx&fpAcc.fig')
saveas(1,'./plt/fx&fpAcc.jpg')


figure(2)
for i = 1:3
    subplot(1,3,i)
    plot(5:0.5:20.5,smooth(Sens(i,:)));hold on;grid on;
    plot(5:0.5:20.5,smooth(fpSens(i,:)));hold off
    title(['data ', num2str(i),' Sensitive']);
    legend('Fixed Point Result','Float Point Result')
    xlabel('SNR')
    ylim([0,1])
end
saveas(2,'./plt/fx&fpSens.fig')
saveas(2,'./plt/fx&fpSens.jpg')
figure(3)
for i = 1:3
    subplot(1,3,i)
    plot(5:0.5:20.5,smooth(FDR(i,:)));hold on;grid on;
    plot(5:0.5:20.5,smooth(fpFDR(i,:)));hold off
    title(['data ', num2str(i),' FDR']);
    legend('Fixed Point Result','Float Point Result')
    xlabel('SNR')
    ylim([0,1])
end
saveas(3,'./plt/fx&fpFDR.fig')
saveas(3,'./plt/fx&fpFDR.jpg')

