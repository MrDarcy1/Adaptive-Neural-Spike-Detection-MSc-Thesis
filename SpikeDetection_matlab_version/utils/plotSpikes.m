function plotSpikes(locations,spikes)
fs = 24414;
plot(0:1/fs:(length(spikes)-1)/fs,spikes);hold on
xlim([0,(length(spikes)-1)/fs])
scatter((locations-1)/fs,spikes(locations),'.');hold off
end

