function [] = fun_plt_cnfdnc_intvl(mu_estimation, sigma_estimation, tit, dim, iter_plot)



m = size(mu_estimation);

cnfdnc_intvl = zeros(2*m(1), m(2));
% 1st row: up confidence interval of 1d
% 2nd row: down confidence interval of 1d
% 3rd row: up confidence interval of 2d
% 4th row: down confidence interval of 2d

for i=1: m(1) % number of dimensiona (2)
    
    for j=1: m(2) % number of points (127)
        mu_temp_dim = mu_estimation(i, j);
        sigma_temp_dim = squeeze(sigma_estimation(j, :, :));
        std_temp = sqrt(sigma_temp_dim(i, i));
        
        up_lim_temp = mu_temp_dim + 2*std_temp;
        down_lim_temp = mu_temp_dim - 2*std_temp;
        
        cnfdnc_intvl( (i-1)*dim + 1, j ) = up_lim_temp;
        cnfdnc_intvl( (i-1)*dim + 2, j ) = down_lim_temp;
    end
end

x_ind = 1: length(cnfdnc_intvl);

switch dim 
    case 1 
        cnfdnc_intvl_d1 = cnfdnc_intvl(1:2, :);
        
        % call plot function 
        funSub_plt_cnfdnt(x_ind, cnfdnc_intvl_d1, mu_estimation, x_sim, tit, 1, iter_plot)
        
    case 2 
        cnfdnc_intvl_d1 = cnfdnc_intvl(1:2, :);
        cnfdnc_intvl_d2 = cnfdnc_intvl(3:4, :);
        
        % call plot function 
        funSub_plt_cnfdnt(x_ind, cnfdnc_intvl_d1, mu_estimation(1,:), cnfdnc_intvl_d2, mu_estimation(2,:), tit, 1, iter_plot, [232, 160, 191]./255, [138, 3, 61]./255, [193, 198, 219]./255, [5, 8, 176]./255)
%         funSub_plt_cnfdnt(x_ind, cnfdnc_intvl_d2, mu_estimation(2,:), tit, 2, iter_plot, [193, 198, 219]./255, [5, 8, 176]./255)
        
% % % % % % % % %         
% % % % % % % % %         figure, 
% % % % % % % % %         hold all 
% % % % % % % % %         plot(x_ind, x_sim(1,:), '--',  'LineWidth', 1, 'Color', [16, 50, 196]./255)
% % % % % % % % %         plot(x_ind, x_sim(2,:), '--', 'LineWidth', 1, 'Color', [201, 24, 10]./255)
% % % % % % % % %         plot(x_ind, mu_estimation(1,:), 'LineWidth', 2, 'Color', [16, 50, 196]./255)
% % % % % % % % %         plot(x_ind, mu_estimation(2,:), 'LineWidth', 2, 'Color', [201, 24, 101]./255)
% % % % % % % % %         legend('x_{d1}', 'x_{d2}', 'estimated /mu_1', 'estimated /mu_2')
% % % % % % % % %         xlim([x_ind(1) x_ind(end)])
% % % % % % % % %         
% % % % % % % % %         str_tit = sprintf('%s -- %d iter', tit, iter_plot);
% % % % % % % % %         title(str_tit)
% % % % % % % % %         hold off 
end


end