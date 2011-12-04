esc = enableSerialControl(esc);
pause(0.1)
esc = setSerialSpeed(esc, 0);
pause(1.5)
esc = setSerialSpeed(esc, 400)
pause(0.05);
[esc dat] = runAdcLog(esc);
esc = setSerialSpeed(esc, 0);

[esc dat] = runAdcLog(esc);

v2 = hex2dec('8000')
v3 = hex2dec('4000')

f2 = dat >= v2;
dat(f2) = dat(f2) - v2;

f3 = dat >= v3;
dat(f3) = dat(f3) - v3;

t = 1:size(dat,2);

subplot(211)
plot(t(f2(1,:)),dat(1,f3(1,:)), ...
    t(f2(2,:)),dat(2,f2(2,:)), ...
    t(f2(3,:)),dat(3,f2(2,:)));
subplot(212)
plot(t(~f2(1,:)),dat(1,~f3(1,:)), ...
    t(~f2(2,:)),dat(2,~f2(2,:)), ...
    t(~f2(3,:)),dat(3,~f2(2,:)));

high = [dat(1,~f3(1,:)); ...
    dat(2,~f2(2,:)); ...
    dat(3,~f2(2,:))];

low = [dat(1,f3(1,:)); ...
    dat(2,f2(2,:)); ...
    dat(3,f2(2,:))];

h(1) = subplot(211);
plot(t(1:2:end),bsxfun(@minus,high',mean(high,1)'),t,zeros(size(t)),'k')
h(2) = subplot(212);
plot(t(1:2:end),bsxfun(@minus,low',mean(low,1)'),t,zeros(size(t)),'k')
linkaxes(h,'x');