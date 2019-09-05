function [noise_data,noise] = addNoise(data,SNRstart, SNRend)
    data=data(:);
%     noise_base=noise_base(:);
    N=length(data);
%     M=length(noise_base);
%     replic=ceil(N/M);
%     noise=repmat(noise_base,replic,1);
%     noise=noise(1:N);
%     noise_p=mean(noise.^2);
%     noise=noise/sqrt(mean(noise.^2));
    signal_power=mean((data-filter(ones(1,16)/16,1,data)).^2);
    delta = (SNRstart-SNRend)/length(data);
%     noise_power=signal_power./(10.^(SNRstart/10))';
    noise_power=signal_power./(10.^([SNRstart:-delta:SNRend+delta]/10));
    noise=sqrt(noise_power).*randn(1,N);
    noise_data=data+noise';
    noise=noise(:);
    
end

