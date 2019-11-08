function [noise_data,noise,numArrival,toa] = addNoisePossion(data,noise_base,SNR,lambda,cells,fs)
    if SNR==inf
        noise_data=data;
        noise=zeros(size(data));
        numArrival=0;
        toa=0;
    else
        data=data(:);
        % noise_power=mean(mean((noise_base).^2)); 
        % noise_base=downsample(reshape(noise_base(:),65,[]),6); 
        for i=1:size(noise_base,2)
            noise_base(:,i)=noise_base(:,i)/max(noise_base(:,i));
        end
        D=64; %duration of noise base
        T=length(data)/fs;
        N=length(data);

        ToA=cell(cells,1);
        for i = 1:cells % toa follows expotiontial distrib.-> can be upgraded to nonhomogenous poission distrib.
            toa=cumsum(exprnd(1/lambda,round(T*(1+lambda)),1));
            ToA{i}=ceil(toa(toa<=T)*fs);
        end

        noise=0;
        % noise_energy=0;
        numArrival=0;
        toa=[];
        for i =1:cells % noise for each cell
            numArrival=numArrival+length(ToA{i});
            toa_temp=ToA{i};
            toa=[toa;toa_temp];
            noise_sel=noise_base(:,randi(1000,1,length(toa_temp)));
            for j =1:length(toa_temp)
                noise_mat=zeros(N,1);
                if toa_temp(j)+D<=N
                    noise_mat(toa_temp(j):toa_temp(j)+D-1)=noise_sel(:,j);
        %             noise_energy=noise_energy+sum(noise_sel(:,j).^2);
                end
                noise=noise+noise_mat;
            end
        end


        % noise_power=noise_energy/(D*numArrival);
        % noise_nomalised=noise/sqrt(noise_power); 
        % noise_nomalised=noise/sqrt(mean(noise.^2)); % the way to evaluate noise power, along all time or only the spike period


        signal_power=mean((data).^2); 
        SNRstart = 20;
        SNRend = 5;
        delta = (SNRstart-SNRend)/length(data);

%         noise_power=signal_power./(10.^([SNRstart:-delta:SNRend+delta]/10))';
        noise_power=signal_power./(10.^(SNR/10));
% 
        noise=sqrt(noise_power).*noise;
        noise_data=data+noise;
    end
end
