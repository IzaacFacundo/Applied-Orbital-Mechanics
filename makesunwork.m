
jdutc = 2458910;

jdut1 = jdutc2jdut1(jdutc, 0);

Tut1 = ((jdut1 - 2451545.0)/36525);

mel = 280.460 + 36000.771*Tut1;

ma = 357.52772333 + 35999.0534*Tut1;

el = mel + 1.914666471*sind(ma) + 0.019994643*sind(2*ma);

r = 1.000140612 - 0.016708617*cosd(ma) - 0.000139589*cosd(2*ma);

obe = 23.439291 - 0.0130042*Tut1;

rtod = [r*cosd(el)  r*cosd(obe)*sind(el) r*sind(obe)*sind(el)]';

rmod = rtod;

zeta = 2306.2181*Tut1 + 0.30188*(Tut1^2) + 0.017998*(Tut1^3); % in arcseconds
theta = 2004.3109*Tut1 - 0.42665*(Tut1^2) - 0.041833*(Tut1^3);
z = 2306.2181*Tut1 + 1.09468*(Tut1^2) + 0.018203*(Tut1^3);

zeta = zeta*pi/648000;
theta = theta*pi/648000;
z = z*pi/648000;

Qmodgcrf = R3(zeta)*R2(-thetalim)*R3(z);

rgcrf = Qmodgcrf*rmod;
posSun = rgcrf;