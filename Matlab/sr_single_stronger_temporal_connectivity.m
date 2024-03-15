clear;clc

% single demo

%% multi-layer modularity calculation

A = FCM.Matrix;
for i=1:length(A)
    A{i} = full(A{i});
end
%模块化参数设置
gamma = 1;
omega = 0.75;

N=length(A{1});
T=length(A);
B=spalloc(N*T,N*T,(N+T)*N*T);
twomu=0;
for s=1:T
    k=sum(A{s});
    twom=sum(k);
    twomu=twomu+twom;
    indx=[1:N]+(s-1)*N;
    B(indx,indx)=A{s}-gamma*k'*k/twom;
end
twomu=twomu+T*omega*N*(T-1);
all2all = N*[(-T+1):-1,1:(T-1)];
B = B + omega*spdiags(ones(N*T,2*T-2),all2all,N*T,N*T);
[S,Q] = genlouvain(B);
Q = Q/twomu
S = reshape(S,N,T);

%% calculate switching rate
% 对于每个节点来说，切换次数除以窗口数量就是切换率
% 节点切换率类似节点的一个属性，只不过必须用多层动态脑网络才能计算出来
% node_num 节点数量
node_num = size(S,1);
% 窗口数量
window_num = size(S,2);
% 对于每一个节点，统计其切换次数，计算切换频率
for i=1:node_num
    node_modules = S(i,:);
    node_diff = diff(node_modules);
    node_sum = sum(node_diff ~= 0);
    node_switching_rate = node_sum / window_num;
    switching_rates(1,i) = node_switching_rate;
end

save('switching_rates.txt','-ascii','switching_rates');
