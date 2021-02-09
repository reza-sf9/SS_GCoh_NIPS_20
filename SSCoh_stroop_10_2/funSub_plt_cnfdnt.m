function [] = funSub_plt_cnfdnt(x_ind, cnfdnc_intvl_data_1, mu_estimation_1, cnfdnc_intvl_data_2, mu_estimation_2, tit, dim_num, iter_plot, patch_color_1, mu_color_1, patch_color_2, mu_color_2)
lbl_size = 60;



figure('units','normalized','outerposition',[0 0 1 1])
hold all

plot(mu_estimation_1(1,:), 'Color', mu_color_1, 'LineWidth', 4);
plot(mu_estimation_2(1,:), 'Color', mu_color_2, 'LineWidth', 4);

patch([x_ind fliplr(x_ind)], [cnfdnc_intvl_data_1(1,:) fliplr(cnfdnc_intvl_data_1(2,:))], patch_color_1, 'LineStyle','none')
patch([x_ind fliplr(x_ind)], [cnfdnc_intvl_data_2(1,:) fliplr(cnfdnc_intvl_data_2(2,:))], patch_color_2, 'LineStyle','none')



xlim([1 x_ind(end)])
str_tit = sprintf('%s -- %d^{st} dimension -- %d Iter ',  tit, dim_num, iter_plot);
% title(str_tit)
xlabel('Trial', 'FontSize',lbl_size)

set(gca,'FontSize',lbl_size)


lgd = legend(' ${\boldmath{x}}_b$', '${\boldmath{x}}_{ic}$');
lgd.FontSize = lbl_size;
set(lgd,'Interpreter','latex');

set(gca, 'XTIck', [20:40:120])

%Returns handles to the patch and line objects
chi=get(gca, 'Children');
%Reverse the stacking order so that the patch overlays the line
set(gca, 'Children',flipud(chi))

hold off

end