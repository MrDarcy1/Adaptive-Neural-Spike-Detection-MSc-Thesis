function [noise_data,noise] = addNoise(data,SNRstart, dummy)
    data=data(:);
    N=length(data);
    signal_power=mean((data-filter(ones(1,16)/16,1,data)).^2);
    noise_power=signal_power./(10.^(SNRstart/10))';
    noise=sqrt(noise_power).*randn(1,N);
    noise_data=data+noise';
    noise=noise(:);
    
end

