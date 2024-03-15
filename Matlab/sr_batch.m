clear;clc;

%% node flexibility, all subjects, fixed parameters.
% 动态脑网络结果文件夹
fcm_root = 'E:\Kaichao_Wu\Documents\3BNTG_day3\Mutilayer\FCM';

subj_names = dir(fcm_root);
% windows 系统设为3
start_id = 3;

% modularity parameters.
gamma = 1;
omega = 1;

for i=start_id:length(subj_names)
    subj_name = subj_names(i).name;
    dfc = load([fcm_root, filesep, subj_name, filesep,'TV_', subj_name, '_FCM.mat']);
    
    
    % modularity calculation, current subject
    A = dfc.FCM.Matrix;
    for k=1:length(A)
        A{k} = full(A{k});
    end

    N=length(A{1}); % T 是节点数目
    T=length(A); % T 是窗口数
    B=spalloc(N*T,N*T,N*N*T+2*N*T);
    twomu=0;
    for s=1:T
        k=sum(A{s});
        twom=sum(k);
        twomu=twomu+twom;
        indx=[1:N]+(s-1)*N;
        B(indx,indx)=A{s}-gamma*k'*k/twom;
    end
    twomu=twomu+2*omega*N*(T-1);
    B = B + omega*spdiags(ones(N*T,2),[-N,N],N*T,N*T);
    [S,Q] = genlouvain(B);
    Q = Q/twomu;
    S = reshape(S,N,T);
    
    % switching rates calculation, current subject
    node_num = size(S,1);
    window_num = size(S,2);
    for j=1:node_num
        node_modules = S(j,:);
        node_diff = diff(node_modules);
        node_sum = sum(node_diff ~= 0);
        node_switching_rate = node_sum / window_num;
        switching_rates(1,j) = node_switching_rate;
    end
    
    switching_rates_all_subjects(i-start_id+1,:) = switching_rates;

end

%% node flexibility, all subjects, multiple parameters.
clear;clc

fcm_root = 'E:\Kaichao_Wu\Documents\3BNTG_day3\Mutilayer\FCM';

subj_names = dir(fcm_root);
start_id = 4;

% modularity parameters.
gammas = [0.9, 1, 1.1];
omegas = [0.5, 0.75, 1];

for i=start_id:length(subj_names)
    subj_name = subj_names(i).name;
    dfc = load([fcm_root, filesep, subj_name, filesep,'TV_', subj_name, '_FCM.mat']);
    
    for gid = 1:length(gammas)       
        for oid= 1:length(omegas)
            
            gamma = gammas(gid);
            omega = omegas(oid);
            
            % modularity calculation, current subject
            A = dfc.FCM.Matrix;
            for k=1:length(A)
                A{k} = full(A{k});
            end

            N=length(A{1});
            T=length(A);
            B=spalloc(N*T,N*T,N*N*T+2*N*T);
            twomu=0;
            for s=1:T
                k=sum(A{s});
                twom=sum(k);
                twomu=twomu+twom;
                indx=[1:N]+(s-1)*N;
                B(indx,indx)=A{s}-gamma*k'*k/twom;
            end
            twomu=twomu+2*omega*N*(T-1);
            B = B + omega*spdiags(ones(N*T,2),[-N,N],N*T,N*T);
            [S,Q] = genlouvain(B);
            Q = Q/twomu;
            S = reshape(S,N,T);

            % switching rates calculation, current subject
            node_num = size(S,1);
            window_num = size(S,2);
            for j=1:node_num
                node_modules = S(j,:);
                node_diff = diff(node_modules);
                node_sum = sum(node_diff > 0);
                node_switching_rate = node_sum / window_num;
                switching_rates(1,j) = node_switching_rate;
            end

            switching_rates_all_subjects(i-start_id+1,gid,oid,:) = switching_rates;
            
        end
    end
    
end

results_folder = 'E:\Kaichao_Wu\Documents\3BNTG_day3\Mutilayer\node_flexibility';
for gid=1:size(switching_rates_all_subjects,2)
    for oid=1:size(switching_rates_all_subjects,3)
        node_flex = squeeze(switching_rates_all_subjects(:,gid,oid,:));
        save([results_folder, filesep, 'node_flex_',num2str(gid),'_',num2str(oid),'.txt'],'-ascii','node_flex');
    end
end