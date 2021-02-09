function [] = fun_plt_simulated_data(xk_sample, dim)

if dim==2
    figure,
    x_color = [148, 28, 1; 18, 17, 18]./255;
    hold all
    plot(1:length(xk_sample), xk_sample(1,:), 'LineWidth', 3, 'Color', x_color(1, :));
    plot(1:length(xk_sample), xk_sample(2,:), 'LineWidth', 3, 'Color', x_color(2, :))
    xlim([1 length(xk_sample)]);
    title('simulated data')
    xlabel('trial')
    legend('1d', '2d')
    hold off
    
elseif dim==1
   
    figure,
    x_color = [148, 28, 1; 18, 17, 18]./255;
    hold all
    plot(1:length(xk_sample), xk_sample(1,:), 'LineWidth', 1.5, 'Color', x_color(1, :));
    xlim([1 length(xk_sample)]);
    title('simulated data')
    xlabel('trial')
    hold off
end

end
