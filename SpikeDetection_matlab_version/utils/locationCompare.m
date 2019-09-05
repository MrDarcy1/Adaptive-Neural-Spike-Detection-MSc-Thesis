function [FP,FN,TP] = locationCompare(locations,interval,spikes)
FP=[];
FN=[];
TP=[];
current_interval=[interval,(1:size(interval,1))'];
for i = 1:length(locations)
    detected=0;
    for j=1:size(current_interval,1)
        if locations(i)>=current_interval(j,1) && locations(i)<=current_interval(j,2)
            TP=[TP spikes(current_interval(j,3))];
            current_interval=[current_interval(1:j-1,:);current_interval(j+1:end,:)];
            detected=1;
            break
        end
    end
    if ~detected
        FN=[FN locations(i)];
    end     
end
try
FP=spikes(current_interval(:,3));
catch
FP=[];
end
end

