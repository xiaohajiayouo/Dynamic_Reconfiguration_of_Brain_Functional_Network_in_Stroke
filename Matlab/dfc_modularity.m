function [ modules_rep, modularity_rep] = dfc_modularity(subj_FCM,gam,ome)


n_sub = size(subj_FCM,1);
T = size(subj_FCM,2);
N = size(subj_FCM,3);
%% Parameters
gamma = gam;
omega = ome;
n_rep = 100;

%% Empty matrices to store data
modularity_rep = zeros(n_sub, n_rep); % Mean modularity
modules_rep = zeros(n_sub,n_rep, N, T); % Module assigment labels    

%%
for sub = 1 : n_sub 
    A = cell(1, T);
    
    % for i=1:length(A)
    %     A{i} = full(A{i});% 将稀疏矩阵转换为稠密矩阵
    % end


    B=spalloc(N*T,N*T,(N+T)*N*T);
    twomu=0;

    for s=1:T
        A{s} = squeeze(subj_FCM(sub,s, :, :) .* (subj_FCM(sub,s, :, :) > 0))
        k    = sum(A{s});
        twom = sum(k);
        twomu = twomu+twom;
        indx=[1:N]+(s-1)*N;
        B(indx,indx) = A{s}-gamma*k'*k/twom;
    end

    twomu=twomu+2*omega*N*(T-1);
    B = B + omega*spdiags(ones(N*T,2),[-N,N],N*T,N*T);
    %%%%% debug %%%%%%%%%%%% 
    [S,Q] = genlouvain(B);
%     Q模块化程�?；S:模块划分结果（所有窗口所有节点）
    Q = Q/twomu % -0.5 - 1 越接�? 模块化程度越�?
    S = reshape(S,N,T); % 把模块划分结果重排列为节点行，层列的矩阵


     for rep = 1 : n_rep
        clc;
        fprintf('repation = %i\n',rep);
        [S,Q] = genlouvain(B);
        Q = Q / twomu;
        S = reshape(S, N, T);

        modularity_rep(sub,rep) = Q;
        modules_rep(sub,rep, :, :) = S;
     end   
end
end