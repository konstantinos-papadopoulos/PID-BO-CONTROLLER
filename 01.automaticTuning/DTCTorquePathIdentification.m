clc
load toKostas
% actual and reference torque for the system identification
% -------------------------------------------------------------------------
Tref     = T_ref_id;
Tact     = Tel_mem ;
Ts       = 2e-5    ;
timeSpan = 0:Ts:length(Tel_mem);
figureInd = 1;
plot(timeSpan,Tref,'b');grid on;hold on
plot(timeSpan,Tact,'r');grid on;hold on
