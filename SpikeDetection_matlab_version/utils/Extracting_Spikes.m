function [data_all_channels,clusters_GT,spikes_all, N, Q] = Extracting_Spikes(folder_raw_data,channels_used)
 
    % Download data
    [data] = SEV2mat(folder_raw_data); % custom function written by Andy's team, converts .sev to .mat
    Fs = data.xWav.fs; % extract sampling rate
    data_all_channels = data.xWav.data(channels_used,:); % all useful raw data
    
    % Important parameters
    nb_of_samples = length(data_all_channels(1,:));
    Q = 1; % number of channels

    
    %% Wave_clus: Spike sorting

    % Save data in .mat files, required in this format for wave_clus
    for i = 1:Q
        data = data_all_channels(i,1:nb_of_samples);
        save(['my_channel_',num2str(i),'.mat'],'data');
        temp = ['my_channel_',num2str(i),'.mat']; % For polytrode
        text(i,1:length(temp)) = temp;
    end
    clear data temp

    % Write polytrode text file - Not required or useful for microwire arrays
    fileID = fopen('polytrode1.txt', 'w');
    formatSpec = '%1s\r\n%2s\r\n';
    for i = 1:Q
        fprintf(fileID,formatSpec,text(i,:));
    end
    fclose(fileID);

    % Setting sampling frequency as parameter
%     param.sr = Fs;

    % Initialising values used in wave_clus iteration indexing
    last_neuron_in_previous_electrode = 0; length_cluster_channel = 0;
    temp1 = 1; temp2 = 0;

    % This for loop serves to organise the detected cluster and spike shapes
    % across the different electrodes into the same matrices: spikes_all =
    % spike shapes; clusters_GT = all neuronal firing times.
    for i = 1:Q
        Get_spikes(['my_channel_',num2str(i),'.mat']) % Wave_clus
        Do_clustering(['my_channel_',num2str(i),'_spikes.mat']) % Wave_clus

        % Storing the data
        % 'where_temp_waveclus_files_are' needs to be the same as Current Folder
        temp_waveclus_file = [pwd,'\times_my_channel_',num2str(i),'.mat'];
        if exist(temp_waveclus_file) ~= 0  % If 0, it means no spikes were detected (or not enough)
            fprintf('Spikes detected in channel \n')
            
            load(['times_my_channel_',num2str(i)]) % Load wave_clus data for the electrode in question
            cluster_class(:,1) = cluster_class(:,1) + last_neuron_in_previous_electrode; % So each detected neuron gets its own number, i.e. channel 1: 1-10, channel 2: 11-15, etc.

            temp2 = temp2 + length(cluster_class(:,1)); % So the indexes line up in the complete list of neurons
            clusters_GT(temp1:temp2,:) = cluster_class; % Store clusters data
            spikes_all(temp1:temp2,:) = spikes; % Store spikes data
            temp1 = temp2+1;

            last_neuron_in_previous_electrode = max(cluster_class(:,1));
            
        else
            fprintf('Spikes not detected in channel \n')
        end
        
        % Deleting surplus files created by wave_clus
        delete(['data_my_channel_',num2str(i),'.dg_01'])
        delete(['data_my_channel_',num2str(i),'.dg_01.lab'])
        delete(['fig2print_my_channel_',num2str(i),'.png'])
    end

    % Clear remaining extra files
%     folder = 'C:\Users\Oscar\OneDrive - Imperial College London\Imperial\PhD\LFP data\Oscars Code\Andrew Jackson\LFPs to Spikes\';
    for i = 1:Q
        delete(['my_channel_',num2str(i),'_spikes.mat'])
        delete(['my_channel_',num2str(i),'.mat'])
        delete(['times_my_channel_',num2str(i),'.mat'])
    end
        delete('spc_log.txt')
        delete('cluster_results.txt')
    for i = [1 4 7 10]
        delete(['fig2print_my_channel_',num2str(i),'a.png'])
    end
    
    % Compensating for offest in neuron nb if some channels don't have neurons
    clusters_GT = sortrows(clusters_GT);
    offset = clusters_GT(1,1)-1;
    clusters_GT(:,1) = clusters_GT(:,1) - offset;

    N = clusters_GT(end,1); % total number of neurons

    
end