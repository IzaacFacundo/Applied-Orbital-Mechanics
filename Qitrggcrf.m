function Q = Qitrggcrf(jdut1, xparc, yparc, dXarc, dYarc, Xarc, Yarc, sarc, dut1) %change time later
   X = arcsec2deg(Xarc);
   Y = arcsec2deg(Yarc);
   dX = arcsec2deg(dXarc);
   dY = arcsec2deg(dYarc);
   xp = arcsec2deg(xparc);
   yp = arcsec2deg(yparc);
   s = arcsec2deg(sarc);
   
   tera = 2*pi*(0.7790572732640 + 1.00273781191135448*(jdut1 - 2451545.0));
   R = R3(-tera);
   
   X = X+dX;
   Y = Y+dY;
   a = 0.5 + 0.125*(X^2 + Y^2);
   PN = [ 1-a*X*2 -a*X*Y X; -a*X*Y 1-a*Y^2 Y; -X -Y 1-a*(X^2+Y^2)];
   PN = deg2rad(PN)*R3(deg2rad(s));
   
   jdtt = (jdut1 - dut1) + (((69.184)/60)/60)/24;
   Ttt = jdtt - ((2451545)/36525);
   sp = -0.000047*Ttt;
   W = R3(-sp)*R2(deg2rad(xp))*R1(deg2rad(yp));
   
   Q = PN*R*W;
end