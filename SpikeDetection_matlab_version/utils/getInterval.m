function [data,spike_locations] = getInterval(data,spike_locations,start,length)
    data = data(start:start+length-1);
    spike_locations = spike_locations(spike_locations>=start);
    spike_locations=spike_locations(spike_locations<start+length)-start+1;
end

