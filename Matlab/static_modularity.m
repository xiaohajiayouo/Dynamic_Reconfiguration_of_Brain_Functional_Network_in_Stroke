%----------------------------------------------------
% Stroke15 : Static normalized modularity calculation

% Kaichao Wu | RMIT-STU
% Last edited: 09-20-2022
%----------------------------------------------------

%% calculate static FC 
clear;
clc;
% loading the denosied parcellation timeseries
denoising_floder = 'H:\Matlab\Work\Stroke15\Mutilayer\ROI_denoising'
% parcellation = {'Network32','Power','Schaefer','AAL'}
parcellation = {'Network32','AAL'}
out_folder = 'H:\Matlab\Work\Stroke15\Mutilayer\Results\FC'

selected_atlas = 1
altas = parcellation{selected_atlas}

denoised_data_path = [denoising_floder,'\',altas,'\'];
save_fc_path = [out_folder,'\',altas,'\'];

subj_path = dir([denoised_data_path,'Sub_0*.txt']);
n_sub = length(subj_path);
for i = 1:n_sub
    subj_name = subj_path(i).name
    denoised_data = load([denoised_data_path,subj_name]);
%     sub_fc(sub_fc<0) = 0;
    A = corr(denoised_data);
    M = atanh(A - diag(diag(A))); % full the diagnoal with 0 and z-score transfered 
    save([save_fc_path,'\',subj_name],'M','-ascii');
    clear M;
end


%% Loading data
clear;
clc;
% parcellation = {'Network32','Power','Schaefer','AAL'}
parcellation = {'Network32','AAL'}
altas_names = {'network32','power','schaefer','aal'}
top_dir = 'H:\Matlab\Work\Stroke15\Mutilayer\Results\FC'
out_dir = 'H:\Matlab\Work\Stroke15\Mutilayer\Results\Static'

for selected_atlas = 1:4
    altas = parcellation{selected_atlas}
    altas_name = altas_names{selected_atlas}
    data_path = [top_dir,'\',altas,'\']
    save_path = [out_dir,'\',altas,'\']

    subj_path = dir([data_path,'Sub_0*.txt'])
    n_sub = length(subj_path)
    n_roi = length(load([data_path,subj_path(1).name]))
    M = zeros(n_sub,n_roi ,n_roi);
    for i = 1:n_sub
        subj_name = subj_path(i).name
        sub_fc = load([data_path,subj_name]);
%         sub_fc(sub_fc<0) = 0; % only positive connectivity left
        M(i,:,:) = sub_fc;
    end

    %% Params

    n_rep = 100;
    n_null = 100;
    n_rew = 1;
    
    density = 0.04:0.01:0.2;
    n_den = length(density)
    n_edges = n_roi*(n_roi-1)
    val_keep_edge = ceil(density.*.5.*n_edges)
    
    q = zeros(n_sub);
    q_null_mean = zeros(n_sub);
    M1 = [1 1 1 1 2 2 2 3 3 3 3 4 4 4 4 4 4 4 5 5 5 5 6 6 6 6 7 7 7 7 8 8]
    %% using the community_louvain_apriori(cla) systematic of negative 
     % the denstiy scheme
    cla_q = zeros(n_sub,n_den);
    cla_q_normalized = zeros(n_sub);    
    for i = 1:n_sub
        i
        sfc = squeeze(M(i,:,:));
        
        At= tril(sfc,-1)
        tril_ind  = tril(true(size(At)),-1);
        val_fc = At(tril_ind);
        sorted_val_fc = sort(abs(val_fc),'descend')
        thres = sorted_val_fc(val_keep_edge)

        for den = 1: n_den
            sfc1 = sfc;
            thre = thres(den);
            sfc1(abs(sfc)<thre) = 0;
            [~,Qt] = community_louvain(sfc1,1,M1,'negative_sym');
            cla_q(i,den) = Qt;
        end 
    end
    save([save_path,altas_name,'_modularity_across_density.mat'], 'cla_q');
   
    
    %% using the gretna modularity calculation
    gretna_q = zeros(n_sub);
    gretna_q_normalized = zeros(n_sub);
    for i = 1:n_sub
        sfc = squeeze(M(i,:,:));
        [A] = gretna_modularity_weight(sfc, '1', n_rep);
        gretna_q(i) = A.modularity_real;
        gretna_q_normalized(i) = A.modularity_zscore;
    end
    
    
    %% customized Static modularity calculation
    for sub = 1: n_sub
        sub
        Qb = 0;
        for rep = 1: n_rep
            A = squeeze(M(sub,:,:)).* squeeze(M(sub, :, :) > 0);
            [~, Qt] = community_louvain(A, 0.9, []); % MD modularit distribution
             if Qt > Qb 
                 Qb = Qt;
             end
             q(sub) = Qb;    
        end
    end
    q_mean = q(:,1)
    % save('mean_power_static_modularity.mat', 'q_mean');

    %% Randomized modularity calculation
    for sub = 1: n_sub
        q_null = zeros(1, n_null);
        sub
        A = squeeze(M(sub,:,:)).* squeeze(M(sub, :, :) > 0);
        if sum(isnan(A(:))) > 0
            q_null = 0;
        else    
            for null = 1 : n_null
                B = randmio_und(A, n_rew);
                Qb = 0;        
                for rep = 1: n_rep 
                    [~, Qt] = community_louvain(B, 1, []);
                    if Qt > Qb 
                        Qb = Qt;
                    end
                end
                q_null(1, null) = Qb;
            end
        end
       q_null_mean(sub) = mean(q_null);       
    end

    %% Normalized modularity calculation
    q_norm = q ./ q_null_mean;
    q_norm1 = q_norm(:,1)
    save([save_path,altas_name,'no_smooth_normalized_modularity.mat'], 'q_norm1');
    save([save_path,altas_name,'no_smooth_mean_modularity.mat'], 'q_mean');
% save('H:\Matlab\Work\Stroke15\Mutilayer\Results\Power\mean_power_normalized_modularity.mat', 'q_norm2');
end

%% average across edge density 
Mild_ind = [1 2 7 10 13 15]
Severe_ind= [3 4 5 6 8 9 11 12 14]
Control_ind = [15:30]
Mild_mean_q = mean(cla_q(Mild_ind ,:))
Severe_mean_q = mean(cla_q(Severe_ind ,:))
Control_mean_q = mean(cla_q(Control_ind ,:))

avg_density = mean(cla_q,2)



