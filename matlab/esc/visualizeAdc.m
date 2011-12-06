
[esc dat] = runAdcLog(esc);

v2 = hex2dec('8000')
v3 = hex2dec('4000')

f2 = dat >= v2;
dat(f2) = dat(f2) - v2;

f3 = dat >= v3;
dat(f3) = dat(f3) - v3;

scale = 3.3 / 1024 * 12.7 / 2.7;

t = 1:size(dat,2);

subplot(211)
plot(t(f3(1,:)),dat(:,f3(1,:)) * scale);
ylim([0 25]);
title('Off period');

subplot(212)
plot(t(~f3(1,:)),dat(:,~f3(1,:))* scale);
ylim([0 25]);
title('On period');
