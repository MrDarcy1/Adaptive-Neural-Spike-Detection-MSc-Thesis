function     visualisation(data,threshold,FP,FN,TP,motion)
fs = 24414;
time = 0:1/fs:(length(data) - 1)/fs;
    if motion==0
        plot(time, data,'b');hold on ;
        plot(time,threshold,'Color','k','LineWidth',1);
        scatter(FN/fs,data(FN),'g+')
        scatter(TP/fs,data(TP),'co');
        scatter(FP/fs,data(FP),'r*');
        hold off
        title('Detection Result')
        legend('Signal','Threshold','FN','TP','FP')
        xlabel('Time Steps')
        ylabel('Amplitude')
        xlim([0,(length(data) - 1)/fs])
    else
        L=500;
        for i = 6000:19:length(data)
            plot(i:i+L,[data(i:i+L)],'b');hold on;
            plot(i:i+L,[threshold(i:i+L)],'k');
            fn=findInRange(FN,[i,i+L]);
            tp=findInRange(TP,[i,i+L]);
            fp=findInRange(FP,[i,i+L]);
            scatter(fn,data(fn),'g+')
            scatter(tp,data(tp),'co')
            scatter(fp,data(fp),'r*')
            hold off;
            xlim([i,i+L+100])
            ylim([-1e-8,3e-8])
            legend('Signal','Threshold','FN','TP','FP')
            pause(0.000001)
        end
    end
            

