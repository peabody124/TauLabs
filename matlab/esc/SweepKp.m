
if ~exist('esc','var');
    esc = EscSerial;
    esc = openPort(esc,'/dev/tty.usbmodemfd131');
end


Kp = [1 2 4 8 10 15 20 25 30 40 50];

for i = 1:length(Kp)
    % Tweak configuration
    esc.configuration.RisingKp = Kp(i);
    esc.configuration.FallingKp = Kp(i);
    esc.configuration.Ki = 50;
    esc.configuration.MaxDcChange = 20;
    esc.configuration.ILim = 1000;
    esc.configuration.InitialStartupSpeed = 30;
    esc = setConfiguration(esc);
    pause(0.1)
    
    [fall_info(i) rise_info(i)] = StepResponse('esc',esc);
end

subplot(221)
plot(Kp, [fall_info.RiseTime]);
xlabel('Kp'); ylabel('Rise time (ms)');

subplot(222);
plot(Kp, [fall_info.Overshoot]);
xlabel('Kp'); ylabel('Overshoot (rpm)');

subplot(223)
plot(Kp, [rise_info.RiseTime]);
xlabel('Kp'); ylabel('Rise time (ms)');

subplot(224);
plot(Kp, [rise_info.Overshoot]);
xlabel('Kp'); ylabel('Overshoot (rpm)');



