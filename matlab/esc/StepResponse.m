if ~exist('esc','var')
    esc = EscSerial;
end

if ~isOpen(esc)
    esc = openPort(esc,'/dev/tty.usbmodemfa141');
end

step_size = 1000;
step_period = 10;
base_speed = 3000;


% Tweak configuration
esc.configuration.RisingKp = 15;
esc.configuration.FallingKp = 1000;
esc.configuration.Ki = 5;
esc.configuration.MaxDcChange = 200;
esc = setConfiguration(esc);
pause(0.1)

% Make sure it's controllable and armed and done with sounds
esc = enableSerialControl(esc);
pause(0.1)
esc = setSerialSpeed(esc, 0);
pause(1.5)


esc = setSerialSpeed(esc, 400)
pause(1);

for i = 400:20:base_speed
    esc = setSerialSpeed(esc, i);
    pause(0.005);
end


esc = enableLogging(esc);

% flush old data
fread(esc.ser, get(esc.ser,'BytesAvailable'), 'uint8');

clf
subplot(211);
h(1) = line('XData', (1:1),'YData', sin(1:1), 'Color', 'b');
h(2) = line('XData', (1:1),'YData', sin(1:1), 'Color', 'g');
legend('speed [rpm]', 'setpoint [rpm]')

t = [];
rpm = [];
setpoint = [];
esc.packet = [];


for i = 1:300
    [self t_ rpm_ setpoint_] = parseLogging(esc);
    t = [t t_];
    rpm = [rpm rpm_];
    if mean(rpm) < 100
        disp('Failed to start or shut down');
        esc = setSerialSpeed(esc, 0);
        break;
    end
    
    setpoint = [setpoint setpoint_];
    
    
    set(h(1), 'XData', t(max(1,end-4000):end), 'YData', rpm(max(1,end-4000):end))
    set(h(2), 'XData', t(max(1,end-4000):end), 'YData', setpoint(max(1,end-4000):end))
    xlim([t(max(1,end-4000+1)) t(end)])
    drawnow;
    
    if mod(i,step_period) == 0 
        esc = setSerialSpeed(esc,base_speed+step_size*(mod(i,step_period*2)==0) );
    end

end

esc = setSerialSpeed(esc, 0);
esc = disableSerialControl(esc);

step_length = 400;

rise = find(setpoint(1:end-1-step_length) == base_speed & setpoint(2:end-step_length) == (base_speed + step_size));
fall = find(setpoint(1:end-1-step_length) == (base_speed + step_size) & setpoint(2:end-step_length) == base_speed);

rise_response = mean(rpm(bsxfun(@plus,rise,(1:step_length)')),2);
fall_response = mean(rpm(bsxfun(@plus,fall,(1:step_length)')),2);

riseinfo = stepinfo(rise_response-base_speed,t(1:step_length));
fallInfo = stepinfo(fall_response-base_speed-step_size,t(1:step_length));

disp(['Rise time was: ' num2str(riseinfo.RiseTime) ' ms and overshoot was ' num2str(riseinfo.Overshoot / step_size * 100) '%']);
disp(['Fall time was: ' num2str(fallInfo.RiseTime) ' ms and overshoot was ' num2str(fallInfo.Overshoot / step_size * 100) '%']);

subplot(212);
plot(t(1:step_length) - t(1), rise_response, t(1:step_length)- t(1), fall_response);
