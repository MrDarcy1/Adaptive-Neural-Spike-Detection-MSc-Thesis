addpath('./data')
addpath('./functions')
addpath('./utils')
%% Load Data
load('realDataWithLFP_1.mat')
load('spike_location_1.mat')
start=randi(length(data)-N);
[data,spike_location]=getInterval(data,spike_location,start,N);

window_length=64;
demean_data=extractMean(data,window_length);    
spike_location=spike_location(spike_location<length(demean_data));

% preprocessed_data=preprocessing(demean_data,'neo',{2,0});
preprocessed_data=preprocessing(demean_data,'aso',{2,0});

% plotProcessedData(spike_location,data,demean_data,preprocessed_data,2)       
 plotProcessedData(spike_location,data,demean_data,preprocessed_data,1,1,{'test','test','test'}) 


%% Thresholding
Mode='rms';
L=64;
c=1:20;
    

Sens=zeros(length(c),1);
FDR=zeros(length(c),1);
Acc=zeros(length(c),1);
parfor j=1:length(c)
    [spikes_detected,threshold,interval,~]=Thresholding_new(data,c(j),L,Mode);
    % Evaluation 
    [FP,FN,TP]=locationCompare(spike_location,interval,spikes_detected);
    Sens(j)= length(TP)/(length(TP)+length(FN)); % found is correct
    FDR(j) = length(FP)/(length(FP)+length(TP)); % not find
    Acc(j) = length(TP)/(length(FP)+length(TP)+length(FN));

end
% figure(32)
% visualisation(preprocessed_data,threshold,FP,FN,TP,0)
% %% Visualisation
    figure(11)
    plot(c,Sens)
    xlabel('c')
%     ylabel('L')
    title('sens')
    figure(21)
    plot(c,FDR)
    xlabel('c')
%     ylabel('L')
    title('fdr')
    figure(32)
    plot(c,Acc)
    xlabel('c')
%     ylabel('L')
    title('acc')

%%
c_aso=15.5;
L_aso=64;
c_neo=12;
L_neo=64;
c_naive=5.4;
%5.6