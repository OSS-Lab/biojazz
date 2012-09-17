function dydt = Ultrasensitive_odes(t, y, rateconstants)

global event_flags;
global event_times

% state vector to node mapping
G0000R = y(1);
TG00000 = y(2);
G0000R_TG00000i00 = y(3);
TG00001 = y(4);
G0000T = y(5);
G0000T_TG00000i00 = y(6);
G0397 = y(7);
G0397_TG00001i00 = y(8);
G0000R_LG0000i00 = y(9);
G0000R_LG0000_TG00000i00 = y(10);
G0000T_LG0000i00 = y(11);
G0000T_LG0000_TG00000i00 = y(12);
LG0000 = y(13);

% constants
fb00 = rateconstants(1);
bb00 = rateconstants(2);
kp00 = rateconstants(3);
fb01 = rateconstants(4);
fu00 = rateconstants(5);
bu00 = rateconstants(6);
fu01 = rateconstants(7);
fu02 = rateconstants(8);
clamp_sink_LG0000 = rateconstants(9);



% expressions
clamp_source_LG0000 = (+(event_flags(1) && ~event_flags(4))*min((t-event_times(1))/1000, 1)*3.33333333333333*4.0+event_flags(4)*max(1-(t-event_times(4))/1000, 0)*3.33333333333333*4.0+(event_flags(2) && ~event_flags(5))*min((t-event_times(2))/1000, 1)*3.33333333333333*4.0+event_flags(5)*max(1-(t-event_times(5))/1000, 0)*3.33333333333333*4.0+(event_flags(3) && ~event_flags(6))*min((t-event_times(3))/1000, 1)*3.33333333333333*4.0+event_flags(6)*max(1-(t-event_times(6))/1000, 0)*3.33333333333333*4.0);

% differential equations for independent species
dydt(1)= + bb00*G0000R_TG00000i00 + kp00*G0000R_TG00000i00 + bb00*G0000R_LG0000i00 + bu00*G0000T - fb00*G0000R*TG00000 - fb01*G0000R*LG0000 - fu00*G0000R ;
dydt(2)= + bb00*G0000R_TG00000i00 + bb00*G0000T_TG00000i00 + kp00*G0397_TG00001i00 + bb00*G0000R_LG0000_TG00000i00 + bb00*G0000T_LG0000_TG00000i00 - fb00*G0000R*TG00000 - fb01*G0000T*TG00000 - fb00*G0000R_LG0000i00*TG00000 - fb01*G0000T_LG0000i00*TG00000 ;
dydt(3)= + fb00*G0000R*TG00000 + bb00*G0000R_LG0000_TG00000i00 + bu00*G0000T_TG00000i00 - bb00*G0000R_TG00000i00 - kp00*G0000R_TG00000i00 - fb01*LG0000*G0000R_TG00000i00 - fu02*G0000R_TG00000i00 ;
dydt(4)= + kp00*G0000R_TG00000i00 + kp00*G0000T_TG00000i00 + bb00*G0397_TG00001i00 + kp00*G0000R_LG0000_TG00000i00 + kp00*G0000T_LG0000_TG00000i00 - fb00*G0397*TG00001 ;
dydt(5)= + bb00*G0000T_TG00000i00 + kp00*G0000T_TG00000i00 + bb00*G0000T_LG0000i00 + fu00*G0000R - fb01*G0000T*TG00000 - fb00*G0000T*LG0000 - bu00*G0000T ;
dydt(6)= + fb01*G0000T*TG00000 + bb00*G0000T_LG0000_TG00000i00 + fu02*G0000R_TG00000i00 - bb00*G0000T_TG00000i00 - kp00*G0000T_TG00000i00 - fb00*LG0000*G0000T_TG00000i00 - bu00*G0000T_TG00000i00 ;
dydt(7)= + bb00*G0397_TG00001i00 + kp00*G0397_TG00001i00 - fb00*G0397*TG00001 ;
dydt(8)= + fb00*G0397*TG00001 - bb00*G0397_TG00001i00 - kp00*G0397_TG00001i00 ;
dydt(9)= + bb00*G0000R_LG0000_TG00000i00 + kp00*G0000R_LG0000_TG00000i00 + fb01*G0000R*LG0000 + bu00*G0000T_LG0000i00 - fb00*G0000R_LG0000i00*TG00000 - bb00*G0000R_LG0000i00 - fu01*G0000R_LG0000i00 ;
dydt(10)= + fb00*G0000R_LG0000i00*TG00000 + fb01*LG0000*G0000R_TG00000i00 + bu00*G0000T_LG0000_TG00000i00 - bb00*G0000R_LG0000_TG00000i00 - kp00*G0000R_LG0000_TG00000i00 - bb00*G0000R_LG0000_TG00000i00 - fu00*G0000R_LG0000_TG00000i00 ;
dydt(11)= + bb00*G0000T_LG0000_TG00000i00 + kp00*G0000T_LG0000_TG00000i00 + fb00*G0000T*LG0000 + fu01*G0000R_LG0000i00 - fb01*G0000T_LG0000i00*TG00000 - bb00*G0000T_LG0000i00 - bu00*G0000T_LG0000i00 ;
dydt(12)= + fb01*G0000T_LG0000i00*TG00000 + fb00*LG0000*G0000T_TG00000i00 + fu00*G0000R_LG0000_TG00000i00 - bb00*G0000T_LG0000_TG00000i00 - kp00*G0000T_LG0000_TG00000i00 - bb00*G0000T_LG0000_TG00000i00 - bu00*G0000T_LG0000_TG00000i00 ;
dydt(13)= + bb00*G0000R_LG0000i00 + bb00*G0000T_LG0000i00 + bb00*G0000R_LG0000_TG00000i00 + bb00*G0000T_LG0000_TG00000i00 + clamp_source_LG0000 - fb01*G0000R*LG0000 - fb00*G0000T*LG0000 - fb01*LG0000*G0000R_TG00000i00 - fb00*LG0000*G0000T_TG00000i00 - clamp_sink_LG0000*LG0000 ;
dydt = dydt(:);
