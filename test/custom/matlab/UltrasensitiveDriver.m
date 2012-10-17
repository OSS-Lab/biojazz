% initial values (free nodes only)
G0000R = 1.00677529813686;
TG00000 = 1.00677529813686;
G0000R_TG00000i00 = 0;
TG00001 = 0;
G0000T = 0;
G0000T_TG00000i00 = 0;
G0397 = 0.1;
G0397_TG00001i00 = 0;
G0000R_LG0000i00 = 0;
G0000R_LG0000_TG00000i00 = 0;
G0000T_LG0000i00 = 0;
G0000T_LG0000_TG00000i00 = 0;
LG0000 = 0.001;
ivalues = [G0000R TG00000 G0000R_TG00000i00 TG00001 G0000T G0000T_TG00000i00 G0397 G0397_TG00001i00 ...
	G0000R_LG0000i00 G0000R_LG0000_TG00000i00 G0000T_LG0000i00 G0000T_LG0000_TG00000i00 LG0000];

% rate constants
fb00= 31.6227766016838;
bb00= 31.6227766016838;
kp00= 515.952796467086;
fb01= 0.001;
fu00= 0.01;
bu00= 1.00451178020473;
fu01= 316.227766016838;
fu02= 3.16227766016838e-07;
clamp_sink_LG0000= 4.0;
rates= [fb00 bb00 kp00 fb01 fu00 bu00 fu01 fu02 clamp_sink_LG0000];

% time interval
t0= 0;
tf= 20000;

% call solver routine 
global event_times;
global event_flags;
[t, y, intervals]= Ultrasensitive_ode_event(@ode23s, @Ultrasensitive_odes, [t0:1:tf], ivalues, odeset('InitialStep', 1e-8, 'AbsTol', 1e-9, 'RelTol', 1e-3, 'MaxStep', 500.0), [0 0 0 0 0 0 0], [500.0], [], [], rates);

% map free node state vector names
G0000R = y(:,1); TG00000 = y(:,2); G0000R_TG00000i00 = y(:,3); TG00001 = y(:,4); G0000T = y(:,5); G0000T_TG00000i00 = y(:,6); G0397 = y(:,7); G0397_TG00001i00 = y(:,8); G0000R_LG0000i00 = y(:,9); G0000R_LG0000_TG00000i00 = y(:,10); 
G0000T_LG0000i00 = y(:,11); G0000T_LG0000_TG00000i00 = y(:,12); LG0000 = y(:,13); 





% issue done message for calling/wrapper scripts
disp('Facile driver script done');

