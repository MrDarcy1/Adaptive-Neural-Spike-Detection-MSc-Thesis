function  plotResultMethods(Sens,FDR,Acc,SNR,lg1,lg2,idx,fig,save,name)
folder='./plt/';
L=length(idx);
for i=1:3 %different data
    figure(fig)
    sens=reshape(Sens(i,idx,:),[L,32,1]);
    for j=1:L
        sens(j,:)=smooth(sens(j,:));
    end
    subplot(1,3,1)
    plot(SNR,sens');
    ylim([0,1])
    legend(lg1{idx})
    title(lg2{i})
    xlabel('SNR')
    ylabel('Sensitive')
%     if save==1
%         saveas(fig,[folder,lg2{i},'_',name,'_Sens.jpg'])
%         saveas(fig,[folder,lg2{i},'_',name,'_Sens.fig'])
%     end
%     fig=fig+1;
    
%     figure(fig)
    fdr=reshape(FDR(i,idx,:),[L,32,1]);
    for j=1:L
        fdr(j,:)=smooth(fdr(j,:));
    end
    subplot(1,3,2)
    plot(SNR,fdr');
    ylim([0,1])
    legend(lg1{idx})
    title(lg2{i})
    xlabel('SNR')
    ylabel('FDR')
%     if save==1
%         saveas(fig,[folder,lg2{i},'_',name,'_FDR.jpg'])
%         saveas(fig,[folder,lg2{i},'_',name,'_FDR.fig'])
%     end
%     fig=fig+1;
    
%     figure(fig)
    acc=reshape(Acc(i,idx,:),[L,32,1]);
    for j=1:L
        acc(j,:)=smooth(acc(j,:));
    end
    subplot(1,3,3)
    plot(SNR,acc');
    ylim([0,1])
    legend(lg1{idx})
    title(lg2{i})
    xlabel('SNR')
    ylabel('Accuracy')
    if save==1
        saveas(fig,[folder,lg2{i},'_',name,'.jpg'])
%         saveas(fig,[folder,lg2{i},'_',name,'.fig'])
    end
    fig=fig+1;

end


end

