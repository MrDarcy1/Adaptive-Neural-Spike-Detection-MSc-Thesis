function [spikes_detected,threshold,interval] = Thresholding_naive(data,c,L)
    N=length(data);
    spikes_detected=[];
    interval=[];
    before=5;
    after=10;
    % initial threshold
    init_threshold=update(data_for_thr(1:L),3); % here 3 is choosen. the data is import to be valid. could be problematice in future.
    threshold=init_threshold*ones(1,L);
    % rest
    hold = 0; %absolute index upto which does not update the threshold
    for i=L+1:N            
        if i-L+1-before>=1
            threshold=[threshold  update(data(i-L+1-before:i-before),Mode,c,threshold(i-1))];
        else %first few data points
            threshold=[threshold  update(data(1:i-before),Mode,c,threshold(i-1))];
        end

        if threshold(i-1)>data(i)
            %find the max location (relative location)
            if i+after<N
                [~,idx]=max(data(i:i+after));
            else %last few data points
                [~,idx]=max(data(i:N));
            end
            spikes_detected=[spikes_detected i+idx-1]; %record spike location
            interval=[interval;[i,i+after]]; % record spike interval
        end          
    end
    threshold=threshold(1:N);
    threshold=threshold(:);
        
    function out=update(data,c) 
        out=mean(data)+c*std(data);
    end
end

