ti_MO = auto_tune_PID_log_parameters(1,1)   ;
tnx_auto = auto_tune_PID_log_parameters(2,1);
tvx_auto = auto_tune_PID_log_parameters(3,1);

x_MO  = tnx_auto + tvx_auto                 ;
y_MO  = tnx_auto * tvx_auto                 ;

num0Gc = 1                       ;
num1Gc = x_MO                    ;
num2Gc = y_MO                    ;
numGc_MO = [num2Gc num1Gc num0Gc];


den0Gc = 0                       ;
den1Gc = ti_MO                   ;
den2Gc = ti_MO*tp6               ;
denGc_MO = [den2Gc den1Gc den0Gc];
Gc_MO = tf(numGc_MO,denGc_MO)    ; 
Gc =  Gc_MO                      ;

% T_net = 0.1*Tp1;
Gp_NET = tf([1],[0.1*Tp1 1]);
Gp_NET_2 = tf([1],[0.3*Tp1 1]);

time = OutputSignal.time;
y_out = OutputSignal.signals.values(:,1);
y_out_with_NET = OutputSignal.signals.values(:,2);
y_out_with_NET_2 = OutputSignal.signals.values(:,3);
plot(time,y_out,'k');hold on; grid on
plot(time,y_out_with_NET,'r');hold on; grid on
plot(time,y_out_with_NET_2,'b');hold on; grid on