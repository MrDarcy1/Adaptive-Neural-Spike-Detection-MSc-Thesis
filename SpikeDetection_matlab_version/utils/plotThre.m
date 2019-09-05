function   plotThre(dataProcessed,threshold,FP,FN,TP,fig,names,save)
    if save==1
        folder='./plotNew/';
        figure(fig)
        visualisation(dataProcessed,threshold,FP,FN,TP,0)
        xlim([1,5000])
        title(names{1})

        saveas(fig,[folder,names{2},'.jpg'])
        saveas(fig,[folder,names{2},'.fig'])
    end
end

