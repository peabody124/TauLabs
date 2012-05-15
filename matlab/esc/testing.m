% Close any open ports and create a new esc device
if ~isempty(instrfind)
fclose(instrfind)
end
port = ls('/dev/tty.usb*'); port(end) = [];
esc = openPort(EscSerial, port)

esc.configuration = getConfiguration(esc)
esc.configuration
esc.configuration.Mode = 0;
esc = setConfiguration(esc)
esc = enableSerialControl(esc)
esc = setSerialSpeed(esc,000);
esc = setSerialSpeed(esc,1000);
pause(1)
esc = setSerialSpeed(esc,3000);
pause(1)
visualizeAdc
getStatus(esc)
esc = setSerialSpeed(esc,400)
pause(1)
[rpm dc current] = sweepRpm(esc);
plot(double(rpm).^2,current,double(rpm).^2,dc)
xlabel('Rpm^2');

b = robustfit(rpm', dc')
Kff = b(2) * 32178 / 1024 * 32
battery_mv = 11000;
Kv = 32178 * 1000 * 2^5 / (battery_mv * Kff)


result.rpm = rpm;
result.dc = dc;
result.cur = current;
result.high_adc = dat(:,~f3(1,:));
result.low_adc = dat(:,f3(1,:));
result.rpm2_power = regress((double(rpm).^2)', double(current)')
