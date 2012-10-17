% initial values (free nodes only)
G0274R = 1e-06;
G0084R0 = 6.30957344480193e-11;
G0084R0_G0274Ri00 = 0;
G0084R1 = 0;
LG0000 = 1e-12;
G0084R0_LG0000i00 = 0;
G1297 = 1.58489319246111e-05;
G1846 = 0.000251188643150957;
G1297_G1846i00 = 0;
LG_0000_x = 0;
ivalues = [G0274R G0084R0 G0084R0_G0274Ri00 G0084R1 LG0000 G0084R0_LG0000i00 G1297 G1846 ...
	G1297_G1846i00 LG_0000_x];

% rate constants
fb00= 3162.27766016838;
bb00= 0.00630957344480193;
kp00= 8.96150501946605;
fb01= 3.16227766016838;
bb01= 1.58489319246112;
bb02= 10;
clamp_sink_LG_0000_x= 1;
rates= [fb00 bb00 kp00 fb01 bb01 bb02 clamp_sink_LG_0000_x];

% time interval
t0= 0;
tf= 100;

% call solver routine 
[t, y]= ode23s(@Template_odes, [t0:0.1:tf], ivalues, odeset('InitialStep', 1e-8, 'AbsTol', 1e-9, 'RelTol', 1e-3, 'MaxStep', 10.0), rates);

% map free node state vector names
G0274R = y(:,1); G0084R0 = y(:,2); G0084R0_G0274Ri00 = y(:,3); G0084R1 = y(:,4); LG0000 = y(:,5); G0084R0_LG0000i00 = y(:,6); G1297 = y(:,7); G1846 = y(:,8); G1297_G1846i00 = y(:,9); LG_0000_x = y(:,10); 






% issue done message for calling/wrapper scripts
disp('Facile driver script done');

