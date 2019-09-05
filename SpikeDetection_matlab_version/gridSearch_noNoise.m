addpath('./data')
addpath('./functions')
%% Load Data
load('realDataWithLFP_1.mat')
load('spike_location_1.mat')
N=100000;
data=-data(1:N);
spike_location=spike_location(spike_location<=N);
%% Preprocessing
% demean
window_length=64;
demean_data=extractMean(data,window_length);
N=N-window_length/2;
% emphasis
method={'aso','neo'};
par={{1,0},{0}};
preprocessed_data=zeros(length(demean_data),length(method));
for i = 1:length(method)
    preprocessed_data(:,i)=preprocessing(demean_data,method{i},par{i});
end

%% Thresholding
Mode='mean';
L=64;
c_aso=20:0.5:25;
c_neo=14:0.5:19;
Sens=zeros(length(c),length(method));
FDR=zeros(length(c),length(method));
Acc=zeros(length(c),length(method));

for j=1:length(c_aso)
    c_(1)=c_aso(j);
    c_(2)=c_neo(j);
    spikes_detected=cell(1,length(method));
    threshold=zeros(size(preprocessed_data));
    interval=cell(1,length(method));
    parfor i = 1:length(method)
        [spikes_detected{i},threshold(:,i),interval{i},~]=thresholding(abs(preprocessed_data(:,i)),c_(i),L,Mode);
    end
    % Evaluation 
    FP=cell(1,length(method));
    TP=cell(1,length(method));
    FN=cell(1,length(method));


    parfor i = 1:length(method)
        [FP{i},FN{i},TP{i}]=locationCompare(spike_location,interval{i},spikes_detected{i});
        Sens(j,i)= length(TP{i})/(length(TP{i})+length(FN{i})); % found is correct
        FDR(j,i) = length(FP{i})/(length(FP{i})+length(TP{i})); % not find
        Acc(j,i) = length(TP{i})/(length(FP{i})+length(TP{i})+length(FN{i}));
    end
end
%% Visualisation
for i = 1:length(method)
    figure(i*3-2)
    plot(Sens(:,i))
    xlabel('c')
%     ylabel('L')
    title('sens')
    figure(i*3-1)
    plot(FDR(:,i))
    xlabel('c')
%     ylabel('L')
    title('fdr')
    figure(i*3)
    plot(Acc(:,i))
    xlabel('c')
%     ylabel('L')
    title('acc')
end
%%
c_aso=15.5;
L_aso=64;
c_neo=12;
L_neo=64;

