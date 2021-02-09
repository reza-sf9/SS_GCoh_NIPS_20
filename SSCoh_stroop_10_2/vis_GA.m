function [] = vis_GA(Bayes_str, dim, iter_plot)

MU_filter = Bayes_str.MU_filter;
SIGMA_filter = Bayes_str.SIGMA_filter;
MU_smoother = Bayes_str.MU_smoother;
SIGMA_smoother = Bayes_str.SIGMA_smoother;



%% plot seperately 
% fun_plt_cnfdnc_intvl(xk_sample, MU_filter, SIGMA_filter, 'Filter Estimation', dim, iter_plot)
fun_plt_cnfdnc_intvl(MU_smoother, SIGMA_smoother, 'Smoother Estimation', dim, iter_plot)


%% subplot of filter smoother
% fun_subPlot_cnfdnc_intv(xk_sample, MU_filter, SIGMA_filter, MU_smoother, SIGMA_smoother, dim, iter_plot)

%% plot summation 
% % if dim==2
% % str_tit_sum = sprintf('summation of x_1 and x_2 --- %d Iter', iter_plot);
% % % fun_plot_sum_2d(xk_sample,MU_filter, MU_smoother, str_tit_sum)
% % end




end