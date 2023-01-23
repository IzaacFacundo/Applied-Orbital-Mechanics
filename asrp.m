function a = asrp(satPos, sunPos, Cr, Amr)
    
    
    rsatsun = sunPos - satPos;
    a = -1/1000 * (1367/299792458) * Cr * Amr * (rsatsun/norm(rsatsun));
end