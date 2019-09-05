addpath('./data')
addpath('./functions')
addpath('./utils')

%% Load Data
load('realDataWithLFP_1.mat')
load('spike_location_1.mat')
N=100000;
data=round((1e6/3)*data(1:N));
% data = data(1:N);
spike_location=spike_location(spike_location<=N);
plotSpikes(spike_location,data)
% spike_location=spike_location(spike_location>N)-N;
fs=24414;
%% Add Noise
load('noise_base.mat');
SNR=inf;
lambda=20; %spike freq
cells=3;
[noise_data,noise]=addNoise(data,3,-5);
% [noise_data,noise,backgroundActNum,backgroundActLoc] = addNoisePossion(data,noise_base,SNR,lambda,cells,fs);

% figure(1)
% subplot(3,1,1)
% plotSpikes(spike_location,data)
% title('Original Data')
% xlabel('Time Steps')
% ylabel('Amplitude')
% % ylim([-3e-4,6e-4])
% xlim([1,5000])
% subplot(3,1,2)
% plotSpikes(spike_location,noise_data)
% xlim([1,5000])
% % ylim([-3e-4,6e-4])
% title(['Signal with Noise'])
% xlabel('Time Steps')
% ylabel('Amplitude')
% subplot(3,1,3)
% plot(noise)
% xlim([1,5000])
% % ylim([-3e-4,6e-4])
% title(['Background Activity - SNR: ',num2str(SNR),'dB'])
% xlabel('Time Steps')
% ylabel('Amplitude')
%% Preprocessing
% demean
window_length=16;
[demean_data,Mean]=extractMean(noise_data,window_length);
% figure(2)
% subplot(2,1,1);
% plot(0:1/fs:(N-1)/fs,data);hold on
% plot(0:1/fs:(length(demean_data)-1)/fs,Mean/16);hold off
% xlim([0,(length(demean_data)-1)/(10*fs)]);
% legend('Original Signl','Mean')
% title('Original Data and Mean')
% xlabel('Time Steps')
% ylabel('Amplitude')
% subplot(2,1,2)
% plotSpikes(spike_location,demean_data)
% xlim([0,(length(demean_data)-1)/(10*fs)]);
% 
% title('Mean Removed Data')
% xlabel('Time Steps')
% ylabel('Amplitude')
%%

N=N-window_length/2;
% emphasis
method={'aso','neo'};
par={{3,0},{3,0}};
preprocessed_data=zeros(length(demean_data),length(method));
for i = 1:length(method)
    preprocessed_data(:,i)=preprocessing(demean_data,method{i},par{i});
    figure(i+2)
    plotSpikes(spike_location,preprocessed_data(:,i));
    title('ASO Emphasised Data k = 3')
    xlabel('Time/s')
    ylabel('Amplitude')
%     ylim([-1e-9,8e-9])
%     xlim([6000,11000])
end

%% Thresholding
Mode='mean';
c=[23,17];
L=[64,64];
% c_aso_rms=16;
% c_neo_rms=12;
% c_aso_mean=23;
% c_neo_mean=17;
spikes_detected=cell(1,length(method));
threshold=zeros(size(preprocessed_data));
interval=cell(1,length(method));
% data_for_thr=zeros(size(preprocessed_data));

for i = 1:length(method)
%     [spikes_detected{i},threshold(:,i),interval{i},dt2]=Thresholding_new(abs(preprocessed_data(:,i)),c(i),L(i),Mode);
    [spikes_detected{i},threshold(:,i),interval{i},dt2]=thresholding(abs(preprocessed_data(:,i)),c(i),L(i),10000);

end

% Evaluation 
FP=cell(1,length(method));
TP=cell(1,length(method));
FN=cell(1,length(method));
Sens=zeros(1,length(method));
FDR=zeros(1,length(method));
Acc=zeros(1,length(method));

for i = 1:length(method)
    [FP{i},FN{i},TP{i}]=locationCompare(spike_location,interval{i},spikes_detected{i});
    Sens(i)= length(TP{i})/(length(TP{i})+length(FN{i})); % found is correct
    FDR(i) = length(FP{i})/(length(FP{i})+length(TP{i})); % not find
    Acc(i) = length(TP{i})/(length(TP{i})+length(FN{i})+length(FP{i}));
end
%% Visualisation
% for i = 1:length(method)
%     figure(3+i+100)
%     visualisation(zeros(size(dt)),threshold(:,i),FP{i},FN{i},TP{i},0)
% end
% figure
visualisation(preprocessed_data(:,1),threshold(:,i),FP{i},FN{i},TP{i},0)
% for i = 1:length(method)
%     figure(3+i+3)
%     visualisation(preprocessed_data(:,i),threshold(:,i),FP{i},FN{i},TP{i},1)
% end
% figure
% plot(preprocessed_data(:,1));hold on
% plot(zeros(size(dt)),'b')
% 
% plot(dt);hold off
%%
% Visualisation(preprocessed_data(:,1),dt1,dt2,FP{i},FN{i},TP{i})

function Visualisation(data,dt1,dt2,FP,FN,TP)

        plot(data,'b');hold on ;
        plot(dt1);
        plot(dt2);
        scatter(FN,data(FN),'r*')
        scatter(TP,data(TP),'co');
        scatter(FP,data(FP),'g+');
        hold off
        title('Detection Result')
        legend('Signal','standard','improved','FN','TP','FP')
        xlabel('Time Steps')
        ylabel('Amplitude')
  end
