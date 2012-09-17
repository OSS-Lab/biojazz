function dydt = Template_odes(t, y, rateconstants)

% state vector to node mapping
G0274R = y(1);
G0084R0 = y(2);
G0084R0_G0274Ri00 = y(3);
G0084R1 = y(4);
LG0000 = y(5);
G0084R0_LG0000i00 = y(6);
G1297 = y(7);
G1846 = y(8);
G1297_G1846i00 = y(9);
LG_0000_x = y(10);

% constants
fb00 = rateconstants(1);
bb00 = rateconstants(2);
kp00 = rateconstants(3);
fb01 = rateconstants(4);
bb01 = rateconstants(5);
bb02 = rateconstants(6);
clamp_sink_LG_0000_x = rateconstants(7);



% expressions
clamp_source_LG_0000_x = (0.5*0.001*1*(1 + square((t-10)/100*2*pi, 50)) + 0.5*0.001*1*(1 + square((t-15)/100*2*pi, 50)) + 0.5*0.001*1*(1 + square((t-20)/100*2*pi, 50)) + 0.5*0.001*1*(1 + square((t-25)/100*2*pi, 50)) + 0.5*0.001*1*(1 + square((t-30)/100*2*pi, 50)))/5;

% differential equations for independent species
dydt(1)= + bb00*G0084R0_G0274Ri00 + kp00*G0084R0_G0274Ri00 - fb00*G0274R*G0084R0 ;
dydt(2)= + bb00*G0084R0_G0274Ri00 + bb01*G0084R0_LG0000i00 - fb00*G0274R*G0084R0 - fb01*G0084R0*LG0000 ;
dydt(3)= + fb00*G0274R*G0084R0 - bb00*G0084R0_G0274Ri00 - kp00*G0084R0_G0274Ri00 ;
dydt(4)= + kp00*G0084R0_G0274Ri00 ;
dydt(5)= + bb01*G0084R0_LG0000i00 - fb01*G0084R0*LG0000 ;
dydt(6)= + fb01*G0084R0*LG0000 - bb01*G0084R0_LG0000i00 ;
dydt(7)= + bb02*G1297_G1846i00 - fb01*G1297*G1846 ;
dydt(8)= + bb02*G1297_G1846i00 - fb01*G1297*G1846 ;
dydt(9)= + fb01*G1297*G1846 - bb02*G1297_G1846i00 ;
dydt(10)= + clamp_source_LG_0000_x - clamp_sink_LG_0000_x*LG_0000_x ;
dydt = dydt(:);

