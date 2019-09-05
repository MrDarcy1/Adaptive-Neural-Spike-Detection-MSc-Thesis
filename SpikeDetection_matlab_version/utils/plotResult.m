function plotResult(Sens,FDR,Acc,SNR,lg1,lg2,idx,fig,save)

for i=idx %different method
    figure(fig)
    fig=fig+1;
    plot(SNR,reshape(Sens(:,i,:),[4,32,1])')
    legend(lg2)
    title(lg1{i})
    xlabel('SNR')
    ylabel('Sensitive')
    
    figure(fig)
    fig=fig+1;
    plot(SNR,reshape(FDR(:,i,:),[4,32,1])')
    legend(lg2)
    title(lg1{i})
    xlabel('SNR')
    ylabel('FDR')
    
    figure(fig)
    fig=fig+1;
    plot(SNR,reshape(Acc(:,i,:),[4,32,1])')
    legend(lg2)
    title(lg1{i})
    xlabel('SNR')
    ylabel('Accuracy')
    if save==1
        saveas(fig,[lg1{i},'.jpg'])
        saveas(fig,[lg1{i},'.fig'])
    end
end
end

