clear;
clc;

path = 'H:\Matlab\Work\Stroke15\Mutilayer\ROI_denoising\Network32\';
savepath = 'H:\Matlab\Work\Stroke15\Mutilayer\Dynamic\Network32\';

step = 1; % window step
window = 50;% window length
n_rep =100;

n_var = 32;

% #########  determinte the number of windows ########
nobs = 200;

if window == nobs
    slides=1;
else
    slides=floor((nobs-window)/step)+1; % how many windows
end
n_win = slides;
% ########################################

subj_path = dir([path,'Sub_0*.txt']);
Subj_FCM = zeros(length(subj_path),n_win,n_var,n_var);

for i =1:length(subj_path) % the number of subjects
    subj_name = subj_path(i).name
    data = load([path,subj_name]);
    FCM = sliding_window_FC(data,window,step,'pos'); 
    Subj_FCM(i,:,:,:) = FCM;
%     input signal data, te window length, step, return the slidling-window based FC matrix
%     Q_f = 0;
%     for j = 1: n_rep
%         [S_t,Q_t] = dfc_modularity(FC);
%         if Q_t>Q_f;
%            Q_f = Q_t;
%            S_f = S_t;
%         end
%     end
      
%     savename = strsplit(subj_name,'.')
%     save([savepath,savename{1,1},'.mat'],'FC')
end
gamma = 0.9;
omega = 1;
[modules_rep,modularity_rep] =  dfc_modularity(Subj_FCM,gamma,omega)
%     [S,Q] = dynamic_modularity_git(FCM);

save(['H:\Matlab\Work\Stroke15\Mutilayer\DFC\Network32\','Modules_rep_50.mat'],'modules_rep')
save(['H:\Matlab\Work\Stroke15\Mutilayer\DFC\Network32\','Modularity_rep_50.mat'],'modularity_rep')
% save(['H:\Matlab\Work\Stroke15\Mutilayer\DFC\','Subj_FCM_4d.mat'],'Subj_FCM')


%% flexibility
nodeflex_rep_mean = Avg_rep_nodeflex(modules_rep);
save(['H:\Matlab\Work\Stroke15\Mutilayer\DFC\Network32\','node_flexibility_tr_50.mat'],'nodeflex_rep_mean')

%% multiple paprmeter
save_path = 'H:\Matlab\Work\Stroke15\Mutilayer\DFC\Modularity_parameters\'

gammas = [0.9, 1, 1.1];
omegas = [0.5, 0.75, 1];
n_gam_ome = length(gammas)+ length(omegas);

for gid = 1:length(gammas)
    for oid= 1:length(omegas)
        gamma = gammas(gid);
        omega = omegas(oid);
        [modules_rep,modularity_rep] =  dfc_modularity(Subj_FCM,gamma,omega);
        save([save_path,['tr_50_Modules_rep_g',num2str(gamma),'_m',num2str(omega),'.mat']],'modules_rep');
        save([save_path,['tr_50_Modularity_rep',num2str(gamma),'_m',num2str(omega),'.mat']],'modularity_rep');
        nodeflex_rep_mean = Avg_rep_nodeflex(modules_rep);
        save([save_path,['tr_50_nodeflex_rep_mean',num2str(gamma),'_m',num2str(omega),'.mat']],'nodeflex_rep_mean');
        clear modules_rep;
        clear modularity_rep;
        clear nodeflex_rep_mean;
    end
end
   
%% load diffeent pa
save_path = 'H:\Matlab\Work\Stroke15\Mutilayer\DFC\Modularity_parameters\'

n_sub = 30
gammas = [0.9, 1, 1.1];
omegas = [0.5, 0.75, 1];
n_gam_ome = length(gammas)+ length(omegas);
Q_max_pa = zeros(n_gam_ome,n_sub);
k =1;
for gid = 1:length(gammas)
    for oid= 1:length(omegas)
        gamma = gammas(gid);
        omega = omegas(oid);
%         load([save_path,['tr_50_Modules_rep_g',num2str(gamma),'_m',num2str(omega),'.mat']])
        load([save_path,['tr_50_Modularity_rep',num2str(gamma),'_m',num2str(omega),'.mat']])
%         load([save_path,['tr_50_nodeflex_rep_mean',num2str(gamma),'_m',num2str(omega),'.mat']])
        for i = 1:30
            Q_max_pa(k,i) = max(modularity_rep(i,:))
        end 
        k = k +1;
    end
end

%% optimal gamma ad omega

anova2(Q_max_qa)
% ttest
for i =1 :9
    [h,p] = ttest2(Q_max_pa(i,1:15)',Q_max_pa(i,15:30)')
end

opt_gamma = gammas(1);
opt_omega = omegas(3);

load([save_path,['nodeflex_rep_mean',num2str(opt_gamma),'_m',num2str(opt_omega),'.mat']])      
        

%% modularity interaction 

% load optimal modualrity and modules assignment
clear;
clc;
top_dir = 'H:\Matlab\Work\Stroke15\Mutilayer\DFC\Modularity_parameters\'
gammas = [0.9, 1, 1.1];
omegas = [0.5, 0.75, 1];
opt_gamma = gammas(1);
opt_omega = omegas(3);
n_gam_ome = length(gammas)+ length(omegas);
n_sub = 30

load([top_dir,['tr_50_Modularity_rep',num2str(opt_gamma),'_m',num2str(opt_omega),'.mat']]) ; 
load([top_dir,['tr_50_Modules_rep_g',num2str(opt_gamma),'_m',num2str(opt_omega),'.mat']]);         
modules_size = size(modules_rep)
n_win = modules_size(4)
Q_max_under_op =  zeros(n_sub);
module_w_interaction = zeros(n_sub,n_win);
module_b_interaction = zeros(n_sub,n_win);
avg_module_w_interaction = zeros(n_sub,n_win);
avg_module_b_interaction = zeros(n_sub,n_win);
n_modules_across_win = zeros(n_sub,n_win);
n_within_modules_across_win = zeros(n_sub,n_win);
n_between_modules_across_win = zeros(n_sub,n_win);
for i = 1:30
     [argvalue,argmax] = max(modularity_rep(i,:));
     max_module = squeeze(modules_rep(i,argmax,:,:));
     
     for win = 1:n_win
         matrix = squeeze(Subj_FCM(i,win,:,:)); 
         modules = max_module(:,win);
         n_modules_across_win(i,win) = length(unique(modules)); % the number of modules each window
         I  = ModularInteraction(matrix, 2, modules);
         [within_module,between_module,avg_within_module,avg_between_module,n_within_modules,n_between_modules] = calculating_avg_w_and_b_interaction(I);
         module_w_interaction(i,win) = within_module;
         module_b_interaction(i,win) = between_module;
         avg_module_w_interaction(i,win) = avg_within_module;
         avg_module_b_interaction(i,win) = avg_between_module;
         n_within_modules_across_win(i,win) = n_within_modules;
         n_between_modules_across_win(i,win) = n_between_modules;
     end
end 
ptients_avg_module_w_interaction = mean( avg_module_w_interaction,2);
ptients_avg_module_b_interaction = mean( avg_module_b_interaction,2);
m_vector = [1,1,0,0,0,0,1,0,0,1,0,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
s_vector = [0,0,1,1,1,1,0,1,1,0,1,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
c_vector = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
m_avg_w_interaction = mean( avg_module_w_interaction(logical(m_vector'),:))
s_avg_w_interaction = mean( avg_module_w_interaction(logical(s_vector'),:))  
c_avg_w_interaction = mean( avg_module_w_interaction(16:30,:))
all_avg_w_interaction = [m_avg_w_interaction',s_avg_w_interaction',c_avg_w_interaction'] 

m_avg_b_interaction = mean( avg_module_b_interaction(logical(m_vector'),:))
s_avg_b_interaction = mean( avg_module_b_interaction(logical(s_vector'),:))  
c_avg_b_interaction = mean( avg_module_b_interaction(16:30,:))
all_avg_b_interaction = [m_avg_b_interaction',s_avg_b_interaction',c_avg_b_interaction'] 

segregtation = (mean(avg_module_w_interaction,2)-mean(avg_module_b_interaction,2))./mean(avg_module_w_interaction,2);

% m_aveg_w_interaction = mean(module_aveg_w_interaction(logical(m_vector'),:))
% s_aveg_w_interaction = mean(module_aveg_w_interaction(logical(s_vector'),:))  
% c_aveg_w_interaction = mean(module_aveg_w_interaction(16:30,:))
% all_aveg_w_interaction = [m_aveg_w_interaction',s_aveg_w_interaction',c_aveg_w_interaction']   
% 
% m_aveg_b_interaction = mean(module_aveg_b_interaction(logical(m_vector'),:))
% s_aveg_b_interaction = mean(module_aveg_b_interaction(logical(s_vector'),:))  
% c_aveg_b_interaction = mean(module_aveg_b_interaction(16:30,:))
% all_aveg_b_interaction = [m_aveg_b_interaction',s_aveg_b_interaction',c_aveg_b_interaction']  
% 
% 
% 
% m_w_interaction = module_aveg_w_interaction(logical(m_vector'),:)
% s_w_interaction = module_aveg_w_interaction(logical(s_vector'),:) 
% c_w_interaction = module_aveg_w_interaction(16:30,:)
% 
% 
% m_b_interaction = module_aveg_b_interaction(logical(m_vector'),:)
% s_b_interaction = module_aveg_b_interaction(logical(s_vector'),:)
% c_b_interaction = module_aveg_b_interaction(16:30,:)
% 
% segregtation2 = (mean(module_aveg_w_interaction,2)-mean(module_aveg_b_interaction,2))./mean(module_aveg_w_interaction,2)
% segregtation = [m_segregation',s_segregation',c_segregation']'

%% the number of modules 
number_of_module = zeros(length(n_sub));
for i  =1:n_sub
    
    c1 = modules_rep(i,:,:,:);
    number_of_module(i)= length(unique(c1));
end

avg_modularity = mean(modularity_rep,2)


num_of_m = [7,7,8,8,8,8,4,6,9,6,6,7,7,6,7,6,10,8,10,9,8,9,8,7,6,8,7,8,9,6]
group = { 'Mild','Mild','Severe','Severe','Severe','Severe','Mild','Severe','Severe','Mild','Severe','Severe','Mild','Severe','Mild',...
    'Control','Control','Control','Control','Control','Control','Control','Control','Control','Control',...
    'Control','Control','Control','Control','Control'};


