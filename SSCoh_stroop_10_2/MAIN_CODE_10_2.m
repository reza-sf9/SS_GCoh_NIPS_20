clc
clear
% clearvars -except ssv count
close all

user_num = 3;

dim = 2;
ch_sim = 5;

new_val = 0;
new_W = new_val;
new_SIM_data = new_val;

coef_length_sim = 10;

coef_noise_W = .1;

Iter = 10;

col_W_1 = 2:5;
col_W_2 = 2:5;
col_W_3 = 2:5;

QQ_sim = [3 0; 0 3].*10^-2;
%% PARAMETER UPDATE
ParamaUpdate.sigma_update = 0; % 1 update, 0 don't update
ParamaUpdate.mu_update = 0;
ParamaUpdate.L_update  = 0;
ParamaUpdate.W_update  = 1;



% load inilial points
fr_vec = [6  40];
sr_vec = [300 256 512 1024];
fr = fr_vec(2); % fre
sr = sr_vec(1); % seg




switch user_num
    case 1
        str_user = '002_jh_a';
    case 2
        str_user = '002_jh_c';
    case 3
        str_user = '011_mm_a';
    case 4
        str_user = '027_ya_a';
end
data_name = str_user;


str_dict = 'C:\Users\REZA_SF\Desktop\NIPS20\nips_code\SSCoh_stroop\data\';

switch ch_sim
    case 5
%         % 2 tasks
        str_name_1 = sprintf('cfg_Init_chNum5_fr6_sr300_%s_winSec2_segNum10_TRIAL.mat',data_name); % 5 channels
        str_name_2 = sprintf('Yk_chNum5_fr6_sr300_%s_winSec2_segNum1_TRIAL',data_name);
        
        % 3 tasks 
%         str_name_1 = sprintf('cfg_Init_TASK_1_2_3_chNum5_fr6_sr300_%s_winSec2_segNum10_TRIAL.mat', data_name); % 5 channels
%         str_name_2 = sprintf('Yk_TASK_1_2_3_chNum5_fr6_sr300_%s_winSec2_segNum1_TRIAL',data_name);
        
    case 20
        str_name_1 = sprintf('cfg_Init_fr%d_sr%d_%s_winSec%d_segNum10_TRIAL.mat', data_name); % 20 channels
        str_name_2 = sprintf('Yk_fr%d_sr%d_%s_winSec%d_segNum1_TRIAL', data_name);
end


str_load_1 = sprintf('%s%s', str_dict, str_name_1);

str_load_2 = sprintf('%s%s', str_dict, str_name_2);


load(str_load_1)
load(str_load_2)


L_eig =  cfg_Init.L;
save('L_eig', 'L_eig')



Yk = Y_k;

y_end = length(Yk);



% call fucntion for calculate initial parameters


if new_W==1
    coef_w = 1;
    w_1 = coef_w*(0 + 1*rand(ch_sim, 1));
    w_1(1,1) = 1;
    w_2 = coef_w*(0 +1*rand(ch_sim, 1));
    w_3 = coef_w*rand(ch_sim, 1);
    W = zeros(ch_sim, 2);
    W(:, 1) = w_1;
    W(:, 2) = w_2;
    W(:, 3) = w_3;
    save('W.txt', 'W', '-ascii')
end

L_init = cfg_Init.L;

% L_r = randn(ch_sim, ch_sim);
% L_i = randn(ch_sim, ch_sim);
% L_init = L_r + L_i*1i;
% 
% L_init = orth(L_init);


% % mu_init = cfg_Init.mu;
incongruent_vec = cfg_Init.incongruent_vec;

K = length(Y_k);         % number of all stepes





if new_SIM_data == 1
    
    %% data simulation
    cnfg_sim.dim = dim;
    cnfg_sim.x0 = [0;0];
    cnfg_sim.F = eye(2);
    cnfg_sim.D = [0;0];
    cnfg_sim.Q = QQ_sim;
    cnfg_sim.incongruent_vec = incongruent_vec;
    cnfg_sim.L = L_init;
    cnfg_sim.W = load('W.txt');
    muu = mean(Y_k).';
    cnfg_sim.mu = muu(1:ch_sim, 1);
    
    % call function to simulate data
    [xk_sample, Yk_sample] = SIM_samp_x_yk(cnfg_sim);
    
    
    str_save = sprintf('X_Y_sim%d.mat', dim);
    save(str_save, 'xk_sample', 'Yk_sample', 'cnfg_sim')
    
    fun_plt_simulated_data(xk_sample, dim)
    
else
%     str_load = sprintf('X_Y_sim%d.mat', dim);
%     load(str_load)
    % % % % % % % % % % % %     fun_plt_simulated_data(xk_sample, dim)
end




%% NEW INITIAL PARAMS
if ParamaUpdate.sigma_update == 1
    init_Param.Q = 100*[1 0; 0 1].*10^-3;
elseif ParamaUpdate.sigma_update == 0 
    init_Param.Q = QQ_sim;
end

if ParamaUpdate.mu_update == 1
    init_Param.mu = cnfg_sim.mu - 4.*(0.5 + 1i*0.5);
elseif ParamaUpdate.mu_update == 0 
    muu = mean(Y_k).';
    init_Param.mu = muu;
end

if ParamaUpdate.L_update == 1
    load('L1');
    init_Param.L  = L1;
elseif ParamaUpdate.L_update == 0
    init_Param.L  = L_init;
end

if ParamaUpdate.W_update == 1
    W_sim = load('W.txt');
    mw = size(W_sim);
    update_coef_W = zeros(mw);
    update_coef_W(col_W_1, 1) =1;
    update_coef_W(col_W_2, 2) =1;
    update_coef_W(col_W_3, 3) =1;

    init_Param.W = W_sim;
    
elseif ParamaUpdate.W_update == 0
    init_Param.W = cnfg_sim.W;
end



%% init param
init_Param.incongruent_vec = incongruent_vec;
init_Param.mu_1_0 = [0; 0];
init_Param.sigma_1_0 = 300*[1 0; 0 1]*10^-3;
init_Param.F = eye(2);
init_Param.D = [0;0];
init_Param.dim = dim;

if dim==1
    init_Param.Q = init_Param.Q(1, 1);
end


PARAM= [];
PARAM{1}=init_Param;
sv_update = [];

step_interval = 1:length(Y_k);
break_con = 0;
thr_q = 10^-5;



for iter=1:Iter
    
    disp(iter)
    disp('---------------------')
    
    param_curr = PARAM{iter};
    
    
    %% filter - smoother
    Bayes = gc_filter_smoother(Yk, param_curr);
    
    Bayes_Step{iter} = Bayes;
    
    %% update parameters
    [updated_param, error_tot] = gc_parameter_update(param_curr, ParamaUpdate, Bayes, Yk, iter, W_sim, update_coef_W);
    
    
   
    
    W_prev = param_curr.W;
    W_update = updated_param.W;
    
    W_diff_prev = W_sim - W_prev;
    W_diff_update = W_sim - W_update;
    
    W_diff_prev_vec = W_diff_prev(:);
    W_diff_update_vec = W_diff_update(:);
    
    vec_diff = W_diff_prev_vec';
    D = sqrt( vec_diff* vec_diff');
    Eug_diff(iter) = D; 
    
    
    
    W_diff_prev_val(iter) = sum(sum(abs(W_diff_prev_vec)));
%     W_diff_update_val(iter) = sum(sum(abs(W_diff_update)));
    
    error_TOT(1, iter) = error_tot;
    
    PARAM{iter+1} = updated_param;
        

    
end

iter_plot = Iter;
Bayes_end = Bayes_Step{iter_plot};
MU_smoother = Bayes_end.MU_smoother;

Param_end = PARAM{iter_plot};
W_end = Param_end.W;

Lambda = zeros(5, length(MU_smoother));
for k=1: length(MU_smoother)
   mu_k =  MU_smoother(:, k);
   e_k = incongruent_vec(k);
   
   x = [mu_k(1) mu_k(2)*e_k 1].';
   
   for ch = 1:5
       W_ch = W_end(ch, :);
       Lambda(ch, k) = exp(W_ch * x);
   end
   
   largest_lamb(1, k) = find(Lambda(:, k) == max(Lambda(:, k)));
   
  
   
    
end

if iter>1
    figure
    plot(Eug_diff), xlabel('iteratoin'), title('Euclidean distance (W_{simulation} and W_{update})')
    xlim([1 length(Eug_diff)])
    
    figure
    plot(W_diff_prev_val), title('compare to real value ')
    xlim([1 length(W_diff_prev_val)])
    
    figure,
    plot(error_TOT), xlabel('iteratoin'), title('error (W_{current} - W_{prev})')
    xlim([1 length(W_diff_prev_val)])
end


plt_ev =1;
if plt_ev==1
yy_vec = largest_lamb;

% marker = {'--or', '--ob', '--ok', '--og', '--oc'};
marker_col = {'r', 'b', 'k', 'g', 'c'};

mrkr_size = 8;

max_val = max(max(Lambda(:)));

figure, hold on 



plot(Lambda(1, :), '-o', 'MarkerFaceColor', marker_col{1}, 'MarkerEdgeColor', marker_col{1}, 'MarkerSize', mrkr_size, 'Color', marker_col{1}); 
plot(Lambda(2, :), '-o', 'MarkerFaceColor', marker_col{2}, 'MarkerEdgeColor', marker_col{2}, 'MarkerSize', mrkr_size, 'Color', marker_col{2} ); 
plot(Lambda(3, :), '-o', 'MarkerFaceColor', marker_col{3}, 'MarkerEdgeColor', marker_col{3}, 'MarkerSize', mrkr_size, 'Color', marker_col{3} ); 
plot(Lambda(4, :), '-o', 'MarkerFaceColor', marker_col{4}, 'MarkerEdgeColor', marker_col{4}, 'MarkerSize', mrkr_size, 'Color', marker_col{4} ); 
plot(Lambda(5, :), '-o', 'MarkerFaceColor', marker_col{5}, 'MarkerEdgeColor', marker_col{5}, 'MarkerSize', mrkr_size, 'Color', marker_col{5} ); 

% legend('\lambda_1','\lambda_2','\lambda_3','\lambda_4','\lambda_5', 'Orientation','horizontal');
% lgd.FontSize = 60;
hold off

yyaxis left
for i=1: length(yy_vec)
    hold on
    plot(i, max_val+(max_val/10), 'o', 'MarkerFaceColor', marker_col{yy_vec(i)}, 'MarkerEdgeColor', marker_col{yy_vec(i)}, 'MarkerSize', mrkr_size ); 
end

hold off






ylabel('\lambda_k','Color', 'k')
xlabel('Trial')
xlim([1 length(Lambda)])
ylim([0 max_val+(max_val/10)])

set(gca, 'XTIck', [0:50:190])

set(gca,'FontSize', 43)
set(gca,'YColor', [0, 0, 0]); % change the color of ytick 

end

if iter>1
    figure
    plot(Eug_diff), xlabel('iteratoin'), title('Euclidean distance (W_{simulation} and W_{update})')
    xlim([1 length(Eug_diff)])
    
    figure
    plot(W_diff_prev_val), title('compare to real value ')
    xlim([1 length(W_diff_prev_val)])
    
    figure,
    plot(error_TOT), xlabel('iteratoin'), title('error (W_{current} - W_{prev})')
    xlim([1 length(W_diff_prev_val)])
end

iter_plot = 1;
param_first = PARAM{iter_plot};
vis_GA(Bayes_Step{iter_plot}, dim, iter_plot);

iter_plot = length(Bayes_Step);
param_end = PARAM{iter_plot};
vis_GA(Bayes_Step{iter_plot}, dim, iter_plot);


mu_first = param_first.mu;
W_first = param_first.W;

mu_end = param_end.mu;
W_end = param_end.W;

for i=1:3
y_label = sprintf('W_%d', i);
w_initial = W_first(:, i);
w_final = W_end(:, i);
vs_vertical_bar_plot(w_initial, w_final, y_label, 60);
end
j


