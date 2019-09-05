function [spikes_detected,threshold,interval,data_for_thr] = thresholding(data,c,L,update_freq)
    data=abs(data);
    N=length(data);
    spikes_detected=[];
    interval=[];
    before=5;
    after=10;
    % initial threshold
    data_for_thr=data;
    init_threshold=update(data_for_thr(1:L),3); % here 3 is choosen. the data is import to be valid. could be problematice in future.
    threshold=init_threshold*ones(1,L);
    %subthreshold=init_threshold*ones(1,L); %not used
    % rest
    hold = 0; %absolute index upto which does not update the threshold
    count=0;
    indicator=0;
    for i=L+1:N 
        if i>hold %not hold
            if threshold(i-1)>data(i) % not detected
                if threshold(i-1)*0.5<=data(i) % avoid the influence before and after a spike.
                    data_for_thr(i-before:i+after)=threshold(i-1)/(c); % notice! only valid when thr=c*metric.
                end

                if count==update_freq
                    if i-L+1-before>=1
                        threshold=[threshold  update(data_for_thr(i-L+1-before:i-before),c)];
                    else %first few data points
                        threshold=[threshold  update(data_for_thr(1:i-before),c)];
                    end
                    count=0;
                else
                    threshold=[threshold threshold(i-1)];
                    count=count+1;
                end
            else %detected
                %find the max location (relative location)
                if i+after<N
                    [~,idx]=max(data(i:i+after));
                else %last few data points
                    [~,idx]=max(data(i:N));
                end
                spikes_detected=[spikes_detected i+idx-1]; %record spike location
                interval=[interval;[i,i+after]]; % record spike interval
                data_for_thr(i-before:i+after)=threshold(i-1)/c; % notice! only valid when thr=c*metric. hold the data for calculate threshold 
                threshold=[threshold threshold(i-1)*ones(1,after+1)]; % hold the threshold
                count=update_freq;
                hold=i+after;
            end          
        end
    end
    threshold=threshold(1:N);
    threshold=threshold(:);

    function out=update(data,c)
        out=c*mean(abs(data));  
    end
end



