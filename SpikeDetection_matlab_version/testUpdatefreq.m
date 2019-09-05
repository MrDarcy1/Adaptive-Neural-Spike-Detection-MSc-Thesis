addpath('./data')
addpath('./functions')
addpath('./utils')

%data setting
N=100000;
fs=24414;
numDataSets=3;
numCombinations=12;

%noise setting
SNR=[10,15,20,inf];
load('noise_base.mat')
lambda=20;
cells=3;

%preprocessing setting
demean_filter_len=16;

%thresholding setting
L=64;
c_aso_mean=25;
update_freq=0:100:10000;

par={'demean','aso','improved mean'};
fig=1;
trial=12;
SenstoPlot=zeros(numDataSets,numCombinations,length(SNR));
FDRtoPlot=zeros(numDataSets,numCombinations,length(SNR));
AcctoPlot=zeros(numDataSets,numCombinations,length(SNR));

for i=1:numDataSets %different data set
%     i=4
    load(['realDataWithLFP_',num2str(i),'.mat'])
    load(['spike_location_',num2str(i),'.mat'])        
    for j=1:length(SNR) %different SNR
        for k=update_freq
            parfor l=1:trial
                if l==1
                    save=1;
                else
                    save=0;
                end
                %%%select interval%%
                start=randi(length(data)-N);
                [data_to_process,spike_location_selected]=getInterval(data,spike_location,start,N);

                %%add noise%%
                [noise_data,noise,~,~] =  addNoisePossion(data_to_process,noise_base,SNR(j),lambda,cells,fs);
                x=[1,5000];
                y=[-3e-4,6e-4];

                %%preprocessing%%         
                %extract mean
                switch par{1}
                    case 'demean'
                        demean_data=extractMean(noise_data,demean_filter_len);
                        spike_location_selected=spike_location_selected(spike_location_selected<length(demean_data));
                    case 'raw'
                        demean_data=data_to_process;
                end
                %emphasis
                switch par{2}
                    case 'aso' 
                        preprocessed_data=preprocessing(demean_data,par{2},{2,0}); %multiRes hop 2
                    case 'neo'
                        preprocessed_data=preprocessing(demean_data,par{2},{2,0}); %multiRes hop 2
                    case 'none'
                        preprocessed_data=demean_data;
                end
              
                %%thresholding%%     
                [spikes_detected,threshold,interval,~]=...
                            thresholding(abs(preprocessed_data),c_aso_mean,L,k);
                        
                %%evaluation%
                [FP,FN,TP]=locationCompare(spike_location_selected,interval,spikes_detected);
                Sens(l) = length(TP)/(length(TP)+length(FN)); % found is correct
                FDR(l) = length(FP)/(length(FP)+length(TP)); % not find
                Acc(l) = length(TP)/(length(TP)+length(FN)+length(FP));

            end
            fig=fig+3;
            SenstoPlot(i,j,k/100+1)=mean(Sens);
            FDRtoPlot(i,j,k/100+1)=mean(FDR);
            AcctoPlot(i,j,k/100+1)=mean(Acc);
            disp(['case: data ',num2str(i),' SNR ',num2str(SNR(j)),'dB undate freq ',num2str(k)])
            disp(['Sens: ',num2str(mean(Sens))])
            disp(['FDR: ',num2str(mean(FDR))])
            disp(['Acc: ',num2str(mean(Acc))])
        end
    end
end    
    for j=1:3 %different SNR
        figure(j)
        subplot(1,3,1)
        acc=reshape(SenstoPlot(j,:,:),12,[]);
        for k=1:3
            acc(k,:)=smooth(acc(k,:));
        end
        plot(update_freq,acc');
        ylim([0,1])
        legend('10','15','20','inf')
        title('Sens')
        subplot(1,3,2)
        acc=reshape(FDRtoPlot(j,:,:),12,[]);
        for k=1:3
            acc(k,:)=smooth(acc(k,:));
        end
        
        plot(update_freq,acc');
                ylim([0,1])

        title('FDR')
        legend('10','15','20','inf')


        subplot(1,3,3)
        acc=reshape(AcctoPlot(j,:,:),12,[]);
        for k=1:3
            acc(k,:)=smooth(acc(k,:));
        end
        plot(update_freq,acc');
                ylim([0,1])

        title('Acc')
        legend('10','15','20','inf')

        saveas(j,['./plt/updateFreq_data_',num2str(j),'.jpg'])
    end
        
