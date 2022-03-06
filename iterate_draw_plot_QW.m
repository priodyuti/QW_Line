  
function iterate_draw_plot_QW(x1,y1,fig,figname,theta,phi1,phi2,color)
  
  
  plot(x1, y1,'color',color,'LineWidth',2);%, '-r','LineWidth',2);
    
  title(['\theta = ', num2str(theta), ', \phi_1 = ',num2str(phi1),', \phi_2 = ',num2str(phi2)]);
  xlabel('k_i','fontweight','bold');
  ylabel('p_i','fontweight','bold');
  %h = legend('SF','Fitted','star','Fitted');
  %set(h,'Interpreter','none')
  %hold off
  saveas(fig,figname)
  clear x1 y1
  
end 

