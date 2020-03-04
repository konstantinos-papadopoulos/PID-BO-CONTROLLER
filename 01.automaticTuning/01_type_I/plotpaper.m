clc
time   = plot_data_u(:,1);
u_auto = plot_data_u(:,2);
u_digi = plot_data_u(:,3);

y_auto = plot_data_y(:,2);
y_digi = plot_data_y(:,3);

figure(1)
plot(time,u_auto,'r',time,u_digi,'k')
grid on


figure(2)
plot(time,y_auto,'r',time,y_digi,'k')
grid on
figure(3)