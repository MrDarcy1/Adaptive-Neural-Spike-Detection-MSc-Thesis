function  plotProcessedData(spike_location,data,demean_data,processedData,fig,save,names)    
    if save==1
        folder='./plotNew/';
        figure(fig)
        subplot(3,1,1)
        plotSpikes(spike_location,data)
        title(names{1})
        xlabel('Time Steps')
        ylabel('Amplitude')
        xlim([1,10000])
        subplot(3,1,2)
        plotSpikes(spike_location,demean_data)
        xlim([1,10000])
        ylim([-4e-5,1e-4])

        title(names{2})
        xlabel('Time Steps')
        ylabel('Amplitude')
        subplot(3,1,3)
        plot(processedData)
        xlim([1,10000])

        title(names{3})
        xlabel('Time Steps')
        ylabel('Amplitude')
%         saveas(fig,[folder,names{4},'.jpg'])
%         saveas(fig,[folder,names{4},'.fig'])
    end
end


