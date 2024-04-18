% uav_visualization 
close all;
t = out.tout;

stepRespons = get(out.logsout, 'step').Values.Data;
step_r = stepRespons(:,1);
step_y = stepRespons(:,2);
step_y_hat = stepRespons(:,3);

sineRespons = get(out.logsout, 'sine').Values.Data;
sine_r = sineRespons(:,1);
sine_y = sineRespons(:,2);
sine_y_hat = sineRespons(:,3);

%% plot style
titlefont = 20;
tickfont = 18;
legendfont = 16;

%% step response
figure();
plot(t, step_r, t, step_y, t, step_y_hat,'-','linewidth',1.5); hold on;
xlabel('Times[s]','interpreter','latex','fontsize',tickfont);
ylabel('Magnitude[m]','interpreter','latex','fontsize',tickfont);
title('Step response', 'fontsize',titlefont);
legend("$r$","$y$","$\hat{y}$",'interpreter','latex','fontsize',legendfont,'location','best')
grid on; grid minor;


%% sinewave tracking
figure();
plot(t, sine_r, t, sine_y, t, sine_y_hat,'-','linewidth',1.5); hold on;
xlabel('Times[s]','interpreter','latex','fontsize',tickfont);
ylabel('Magnitude[m]','interpreter','latex','fontsize',tickfont);
title('400 Hz Sinewave tracking', 'fontsize',titlefont);
legend("$r$","$y$","$\hat{y}$",'interpreter','latex','fontsize',legendfont,'location','best')
grid on; grid minor;

% figure();
% colors = colormap(jet(length(t))); 
% for i = 1:100:length(t)
%   plot(x(i), y(i),'.','linewidth',1.5,"color",colors(i,:)); hold on;
% end
% c=colorbar('Ticks',linspace(0,1,11),'TickLabels',linspace(0,t(end),11));
% c.Label.String = 'Time[s]';
% 
% % plot(x, y,'-','linewidth',1.5);
% xlabel('x[m]','interpreter','latex','fontsize',12);
% ylabel('y[m]','interpreter','latex','fontsize',12);
% title('Position[ENU]');
% legend("2D position(XY)",'interpreter','latex','fontsize',12,'location','best')
% grid on; grid minor;
