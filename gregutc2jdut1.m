function jdut1 = gregutc2jdut1(time,dut1)
    mon = time(1);
    day = time(2);
    year = time(3);
    hr = time(4);
    min = time(5);
    sec = time(6);

    jdutc = (367*year) - fix(7*(year+fix((mon+9)/12))/4) + fix((275*mon)/9) + day +1721013.5 + (1/24)*(hr+(1/60)*(min+(sec/60)));
    jdut1 = jdutc + dut1;
end