addpath('./data')
addpath('./functions')
addpath('./utils')

%data setting
N=100000;
fs=24414;
numDataSets=3;
numCombinations=12;

%noise setting
SNR=[5:0.5:20,inf];
load('noise_base.mat')
lambda=20;
cells=3;

%preprocessing setting
demean_filter_len=16;

%thresholding setting
L=64;
c_none_naive=8;
c_neo_naive=55;
c_aso_naive=55;

c_none_median=7;
c_neo_median=55;
c_aso_median=55;

c_none_mean=6;
c_neo_mean=20;
c_aso_mean=25;

c_none_rms=5;
c_neo_rms=15;
c_aso_rms=16;

c_raw_median=6;
c_raw_mean=5;
c_raw_rms=3; % acc about .25 ~.4

c=[c_none_naive,c_neo_naive,c_aso_naive,...
   c_none_median,c_neo_median,c_aso_median,...
   c_none_mean,c_neo_mean,c_aso_mean,...
   c_none_rms,c_neo_rms,c_aso_rms,...
   c_none_mean,c_neo_mean,c_aso_mean,...
   c_none_rms,c_neo_rms,c_aso_rms,...
   c_none_median,c_neo_median,c_aso_median,...
   c_none_mean,c_neo_mean,c_aso_mean,...
   c_none_rms,c_neo_rms,c_aso_rms,...
   c_raw_median,c_raw_mean,c_raw_rms]; %c_naive


parameters={{'demean','none','naive median'},{'demean','neo','naive median'},{'demean','aso','naive median'},...
            {'demean','none','median'},{'demean','neo','median'},{'demean','aso','median'}...
            {'demean','none','mean'},{'demean','neo','mean'},{'demean','aso','mean'},...
            {'demean','none','rms'},{'demean','neo','rms'},{'demean','aso','rms'},...
            {'demean','none','naive mean'},{'demean','neo','naive mean'},{'demean','aso','naive mean'},...
            {'demean','none','naive rms'},{'demean','neo','naive rms'},{'demean','aso','naive rms'},...
            {'demean','none','improved median'},{'demean','neo','improved median'},{'demean','aso','improved median'},...
            {'demean','none','improved mean'},{'demean','neo','improved mean'},{'demean','aso','improved mean'},...
            {'demean','none','improved rms'},{'demean','neo','improved rms'},{'demean','aso','improved rms'},...
            {'raw','none','median'},{'raw','none','mean'},{'raw','none','rms'}};
case_to_test=28:30;
fig=1;
trial=12;
% SenstoPlot=zeros(numDataSets,numCombinations,length(SNR));
% FDRtoPlot=zeros(numDataSets,numCombinations,length(SNR));
% AcctoPlot=zeros(numDataSets,numCombinations,length(SNR));
SNRtoPlot=5:0.5:20.5;
for i=1:numDataSets %different data set
%     i=4
    load(['realDataWithLFP_',num2str(i),'.mat'])
    load(['spike_location_',num2str(i),'.mat'])
    for k=case_to_test %different method
        par=parameters{k};         
        for j=1:length(SNR) %different SNR
%             j=length(SNR);
            for l=1:trial
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
                names={['data ',num2str(i)],['data_',num2str(i),'_SNR_',num2str(SNR(j)),'dB']};
                plotNoise(spike_location_selected,data_to_process,noise_data,noise,SNR(j),x,y,fig,names,save)     

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
                names={['data ',num2str(i)],[par{1}, ' data'],[par{2},' processed data'],...
                       ['data_',num2str(i),'_SNR_',num2str(SNR(j)),'dB_',par{1},'_',par{2}]};
%                  plotProcessedData(spike_location_selected,data_to_process,demean_data,preprocessed_data,fig+1,save,names) 

                
                %%thresholding%%     
                [spikes_detected,threshold,interval,~]=...
                            Thresholding_new(abs(preprocessed_data),c(k),L,par{3});
                        
                %%evaluation%
                [FP,FN,TP]=locationCompare(spike_location_selected,interval,spikes_detected);
                Sens(l) = length(TP)/(length(TP)+length(FN)); % found is correct
                FDR(l) = length(FP)/(length(FP)+length(TP)); % not find
                Acc(l) = length(TP)/(length(TP)+length(FN)+length(FP));
                names={['data ',num2str(i),' SNR ',num2str(SNR(j)),' dB ',par{1},' ',par{2},' ',par{3}],...
                       ['data_',num2str(i),'_SNR_',num2str(SNR(j)),' dB_',par{1},'_',par{2},'_',par{3}]};
                 plotThre(preprocessed_data,threshold,FP,FN,TP,fig+2,names,save)
% 
%                 disp(['case: data ',num2str(i),' SNR ',num2str(SNR(j)),'dB ',par{1},' ',par{2},' ',par{3}])
%                 disp(['Sens: ',num2str(Sens(l))])
%                 disp(['FDR: ',num2str(FDR(l))])
%                 disp(['Acc: ',num2str(Acc(l))])
            end
            fig=fig+3;
            SenstoPlot(i,k,j)=mean(Sens);
            FDRtoPlot(i,k,j)=mean(FDR);
            AcctoPlot(i,k,j)=mean(Acc);
            disp(['case: data ',num2str(i),' SNR ',num2str(SNR(j)),'dB ',par{1},' ',par{2},' ',par{3}])
            disp(['Sens: ',num2str(mean(Sens))])
            disp(['FDR: ',num2str(mean(FDR))])
            disp(['Acc: ',num2str(mean(Acc))])
        end
    end
end

par=parameters{1};
lg1={[par{1},' ',par{2},' ',par{3}]};
for i =2:length(parameters)
    par=parameters{i};
    lg1=[lg1,[par{1},' ',par{2},' ',par{3}]];
end
lg2={'data 1','data 2','data 3'};
for i =2:3
    lg2=[lg2,['data ',num2str(i)]];
end
% 
plotResultMethods(SenstoPlot,FDRtoPlot,AcctoPlot,SNRtoPlot,lg1,lg2,[1,2,3],fig,1,'naive_median')
plotResultMethods(SenstoPlot,FDRtoPlot,AcctoPlot,SNRtoPlot,lg1,lg2,[13,14,15],fig,1,'naive_mean')
plotResultMethods(SenstoPlot,FDRtoPlot,AcctoPlot,SNRtoPlot,lg1,lg2,[16,17,18],fig,1,'naive_rms')

plotResultMethods(SenstoPlot,FDRtoPlot,AcctoPlot,SNRtoPlot,lg1,lg2,[4,5,6],fig,1,'standard_median')
plotResultMethods(SenstoPlot,FDRtoPlot,AcctoPlot,SNRtoPlot,lg1,lg2,[7,8,9],fig,1,'standard_mean')
plotResultMethods(SenstoPlot,FDRtoPlot,AcctoPlot,SNRtoPlot,lg1,lg2,[10,11,12],fig,1,'standard_rms')

plotResultMethods(SenstoPlot,FDRtoPlot,AcctoPlot,SNRtoPlot,lg1,lg2,[19,20,21],fig,1,'improved_median')
plotResultMethods(SenstoPlot,FDRtoPlot,AcctoPlot,SNRtoPlot,lg1,lg2,[22,23,24],fig,1,'improved_ mean')
plotResultMethods(SenstoPlot,FDRtoPlot,AcctoPlot,SNRtoPlot,lg1,lg2,[25,26,27],fig,1,'improved_rms')


% plotResultMethods(SenstoPlot,FDRtoPlot,AcctoPlot,SNRtoPlot,lg1,lg2,[1,13,16,4,7,10],fig,1,'none')%,19,22,25],fig,1,'none')
% plotResultMethods(SenstoPlot,FDRtoPlot,AcctoPlot,SNRtoPlot,lg1,lg2,[2,14,17,5,8,11],fig,1,'neo')%,20,23,26],fig,1,'neo')
% plotResultMethods(SenstoPlot,FDRtoPlot,AcctoPlot,SNRtoPlot,lg1,lg2,[3,15,18,6,9,12],fig,1,'aso')%,21,24,27],fig,1,'aso')
% 
plotResultMethods(SenstoPlot,FDRtoPlot,AcctoPlot,SNRtoPlot,lg1,lg2,[3,6,21],fig,1,'aso_median')
plotResultMethods(SenstoPlot,FDRtoPlot,AcctoPlot,SNRtoPlot,lg1,lg2,[15,9,24],fig,1,'aso_mean')
plotResultMethods(SenstoPlot,FDRtoPlot,AcctoPlot,SNRtoPlot,lg1,lg2,[18,12,27],fig,1,'aso_rms')
% 
plotResultMethods(SenstoPlot,FDRtoPlot,AcctoPlot,SNRtoPlot,lg1,lg2,[3,15,18],fig,1,'naive')
plotResultMethods(SenstoPlot,FDRtoPlot,AcctoPlot,SNRtoPlot,lg1,lg2,[6,9,12],fig,1,'standard')
plotResultMethods(SenstoPlot,FDRtoPlot,AcctoPlot,SNRtoPlot,lg1,lg2,[21,24,27],fig,1,'improved')









% for i=1:4 %different data set
%     figure(fig)
%     fig=fig+1;
%     plot(SNRtoPlot,reshape(SenstoPlot(i,:,:),[9,32,1])')
%     legend(lg1)
%     title(lg2{i})
%     xlabel('SNR')
%     ylabel('Sensitive')
%     
%     figure(fig)
%     fig=fig+1;
%     plot(SNRtoPlot,reshape(FDRtoPlot(i,:,:),[9,32,1])')
%     legend(lg1)
%     title(lg2{i})
%     xlabel('SNR')
%     ylabel('FDR')
%     
%     figure(fig)
%     fig=fig+1;
%     plot(SNRtoPlot,reshape(AcctoPlot(i,:,:),[9,32,1])')
%     legend(lg1)
%     title(lg2{i})
%     xlabel('SNR')
%     ylabel('Accuracy')
% %     saveas(fig,[lg2{i},'.jpg'])
% %     saveas(fig,[lg2{i},'.fig'])
% end
% 
% for i=1:9 %different method
%     figure(fig)
%     fig=fig+1;
%     plot(SNRtoPlot,reshape(SenstoPlot(:,i,:),[4,32,1])')
%     legend(lg2)
%     title(lg1{i})
%     xlabel('SNR')
%     ylabel('Sensitive')
%     
%     figure(fig)
%     fig=fig+1;
%     plot(SNRtoPlot,reshape(FDRtoPlot(:,i,:),[4,32,1])')
%     legend(lg2)
%     title(lg1{i})
%     xlabel('SNR')
%     ylabel('FDR')
%     
%     figure(fig)
%     fig=fig+1;
%     plot(SNRtoPlot,reshape(AcctoPlot(:,i,:),[4,32,1])')
%     legend(lg2)
%     title(lg1{i})
%     xlabel('SNR')
%     ylabel('Accuracy')
% %     saveas(fig,[lg1{i},'.jpg'])
% %     saveas(fig,[lg1{i},'.fig'])
% end

% for l=1:10
%     for i=1:4 %different data set
%         %load data
%         load(['realDataWithLFP_',num2str(i),'.mat'])
%         load(['spike_location_',num2str(i),'.mat'])
%         start=randi(length(data));
%         N=10000;
%         [data,spike_location]=getInterval(data,spike_location,start,N);
%         for j=1:length(SNR) %different SNR
%             %add noise
%             [noise_data,noise,~,~] =  addNoisePossion(data,noise_base,SNR(j),lambda,cells,fs);
%             x=[1,5000];
%             y=[-3e-4,6e-4];
%             names={['data ',num2str(i)],['data_',num2str(i),'_SNR_',num2str(SNR(j)),'dB']};
%             plotNoise(spike_location,data,noise_data,noise,SNR(j),x,y,fig,names,save)     
%             fig=fig+1;
%             for k=1:9 %different method
%                 par=parameters{k};         
%                 %%preprocessing%%         
%                 %extract mean
%                 switch par{1}
%                     case 'demean'
%                         demean_data=extractMean(noise_data,demean_filter_len); 
%                     case 'raw'
%                         demean_data=data;
%                 end
%                 %emphasis
%                 switch par{2}
%                     case 'aso' 
%                         preprocessed_data=preprocessing(demean_data,par{2},{2,0}); %multiRes hop 2
%                     case 'neo'
%                         preprocessed_data=preprocessing(demean_data,par{2},{2,0}); %multiRes hop 2
%                     case 'none'
%                         preprocessed_data=demean_data;
%                 end
%                 save=1;
%                 names={['data ',num2str(i)],[par{1}, ' data'],[par{2},' processed data'],...
%                        ['data_',num2str(i),'_SNR_',num2str(SNR(j)),'dB_',par{1},'_',par{2}]};
%                 plotProcessedData(spike_location,data,demean_data,preprocessed_data,fig,save,names) 
%                 fig=fig+1;
%                 %%thresholding%%     
%                 [spikes_detected,threshold,interval,~]=...
%                             Thresholding_new(abs(preprocessed_data),c(k),L,par{3});
%                 %%evaluation%
%                 [FP,FN,TP]=locationCompare(spike_location,interval,spikes_detected);
%                 Sens = length(TP)/(length(TP)+length(FN)); % found is correct
%                 FDR = length(FP)/(length(FP)+length(TP)); % not find
%                 Acc = length(TP)/(length(TP)+length(FN)+length(FP));
%                 names={['data ',num2str(i),' SNR ',num2str(SNR(j)),'dB ',par{1},' ',par{2},' ',par{3}],...
%                        ['data_',num2str(i),'_SNR_',num2str(SNR(j)),'dB_',par{1},'_',par{2},'_',par{3}]};
%                 plotThre(preprocessed_data,threshold,FP,FN,TP,fig,names,save)
%                 fig=fig+1;
%                 disp(['case: data ',num2str(i),' SNR ',num2str(SNR(j)),'dB ',par{1},' ',par{2},' ',par{3}])
%                 disp(['Sens: ',num2str(Sens)])
%                 disp(['FDR: ',num2str(FDR)])
%                 disp(['Acc: ',num2str(Acc)])
% 
%             end
%         end
%     end
%     
% end
