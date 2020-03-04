figure(1)
plot1=plot(Time,data.signals.values(:,1),Time,data.signals.values(:,3),Time,data.signals.values(:,4),Time,data.signals.values(:,2),Time,data.signals.values(:,5),Time,data.signals.values(:,6),'LineWidth',1);
set(plot1(1),'Color','black','LineStyle','--','DisplayName','I_ref');
set(plot1(4),'Color','black','DisplayName','I_{FEED}');
set(plot1(2),'DisplayName','I_{FEED} +\alpha');
set(plot1(3),'DisplayName','I_{FEED} -\alpha');
set(plot1(5),'Color','Black','LineStyle','-.','DisplayName','V_{NET}');
set(plot1(6),'Color','Black','LineStyle','--','DisplayName','I_{LOAD}');
h1 = legend('show');
set(h1,'Location','NorthWest','FontSize',10,'FontName','Arial','FontName','Arial');
xlabel('Time [s]')
ylabel('Control values [p.u.]')

figure(2)
plot2=plot(Time,data.signals.values(:,1),Time,data.signals.values(:,3),Time,data.signals.values(:,4),Time,data.signals.values(:,2),'LineWidth',1);
set(plot2(1),'Color','black','LineStyle','--','DisplayName','I_ref');
set(plot2(4),'Color','black','DisplayName','I_{FEED}');
set(plot2(2),'DisplayName','I_{FEED} +\alpha');
set(plot2(3),'DisplayName','I_{FEED} -\alpha');
h1 = legend('show');
set(h1,'Location','NorthEast','FontSize',10,'FontName','Arial','FontName','Arial');
xlabel('Time [s]')
ylabel('Control values [p.u.]')
xlim([0.3 0.7])
ylim([0.9 1.2])