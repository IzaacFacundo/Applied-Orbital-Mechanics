frac = [78517110 4026510 4026510 124822050 132875100 122808810 161877990 3189330 116948430 174575490 139413210 191179890];
pfrac = [.0833 .25 .25 .25 .25 .25 .083 .25 .083 .083 .25 .083];



s = 0:100000:200000000;
so = 100000000;
m = 5;
P5 = exp(-(s./so).^m);
m = 10;
P10 = exp(-(s./so).^m);
m = 20;
P20 = exp(-(s./so).^m);



fig = 1;
figure(fig);
plot(s,P5, 'Linewidth', 3)
hold on
plot(s,P10, 'Linewidth', 3)
plot(s,P20, 'Linewidth', 3)
scatter(frac,pfrac)
xlabel('Stress/Median Stress (Pa)')
ylabel('Probability')
legend('m = 5', 'm=10', 'm=20')
title('Weibull Distributions')
