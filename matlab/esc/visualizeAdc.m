
[esc dat] = runAdcLog(esc);

v2 = hex2dec('8000');
v3 = hex2dec('4000');

f2 = dat >= v2;
dat(f2) = dat(f2) - v2;

f3 = dat >= v3;
dat(f3) = dat(f3) - v3;

scale = 3.3 / 4096 * 12.7 / 2.7;

status = getStatus(esc);
mean(max(dat(2:end,:),[],1) * scale) / (double(status.DutyCycle) / 100)

pwm_rate = 80000;
t = (1:size(dat,2)) / pwm_rate * 1000;

f3(1,:) = [];

subplot(211)
plot(t(f3(1,:)),dat(:,f3(1,:)) * scale);
title('Off period');
ylim([0 15]);

subplot(212)
plot(t(~f3(1,:)),dat(:,~f3(1,:))* scale);
ylim([0 15]);
title('On period');
xlabel('Time (ms)')