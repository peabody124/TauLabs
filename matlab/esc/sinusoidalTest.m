function [t rpm setpoint] = sinusoidalTest(varargin)

params.esc = [];
params.port = '/dev/tty.usbmodemfa141';
params.dt_ms = 1;
params.base_speed_rpm = 4000;
params.sin_amp_rpm = 50;
params.start_freq_hz = 0.1;
params.end_freq_hz = 20;
params.test_time_s = 100;
params = parseVarArgs(params,varargin{:});

% If no esc object passed in, create one
esc = params.esc;
if isempty(params.esc)
    esc = EscSerial;
end

% Open port if not open
if ~isOpen(esc)
    esc = openPort(esc,'/dev/tty.usbmodemfa141');
    opened = true;
else
    opened = false;
end

% Make sure it's controllable and armed and done with sounds
esc = enableSerialControl(esc);
pause(0.1)
esc = setSerialSpeed(esc, 0);
pause(3.5)

% Start up ESC then ramp it up
esc = setSerialSpeed(esc, 550);
pause(4);
for i = 550:5:params.base_speed_rpm
    esc = setSerialSpeed(esc, i);
    pause(0.01);
end

% flush old data and enable logging 
if get(esc.ser,'BytesAvailable') > 0
    fread(esc.ser, get(esc.ser,'BytesAvailable'), 'uint8');
end
esc = enableLogging(esc);
esc.packet = [];

t = [];
rpm = [];
setpoint = [];

low = 0;

tic;
a = [];
while toc < params.test_time_s
    pause(params.dt_ms / 1000);
    freq = params.start_freq_hz + toc * (params.end_freq_hz - params.start_freq_hz) / params.test_time_s;
    speed = sin(2*pi*freq*toc);
    speed = params.base_speed_rpm + params.sin_amp_rpm * speed;
    esc = setSerialSpeed(esc, speed);
    
    [esc t_ rpm_ setpoint_] = parseLogging(esc);    
    t = [t t_];
    rpm = [rpm rpm_];
    if mean(rpm) < 100
        disp('Failed to start or shut down');
        esc = setSerialSpeed(esc, 0);
        break;
    end
    setpoint = [setpoint setpoint_];
end

% Shut down controller and resume PWM control.
esc = setSerialSpeed(esc, 0);
esc = disableSerialControl(esc);
esc = disableLogging(esc);
% If we opened the port close it again
if opened
    esc = closePort(esc);
end

clf
subplot(211);
plot((t-t(1))/1000,setpoint,(t-t(1))/1000,rpm + 2* params.sin_amp_rpm);
legend('setpoint [rpm]','speed [rpm]')
title('Input and Output (intentionally offset)')
ylim([params.base_speed_rpm-params.sin_amp_rpm*2  params.base_speed_rpm + 4* params.sin_amp_rpm]);
xlabel('Time (s)')
ylabel('RPM');

F = 1000/mean(diff(t));
rpm = medfilt1(rpm,3);
setpoint = medfilt1(setpoint,3);

for i = 1:floor(length(rpm)/1000)
    input = fft(setpoint(1+(i-1)*1000:1000+(i-1)*1000)-mean(setpoint(1+(i-1)*1000:1000+(i-1)*1000)));
    output = fft(rpm(1+(i-1)*1000:1000+(i-1)*1000)-mean(rpm(1+(i-1)*1000:1000+(i-1)*1000)));
    [~,idx] = max(abs(input));
    freq(i) = idx; % Because sample length is 1000 and sampling is 1k index is freq
    phase_shift(i) = angle(output(idx)/input(idx));
    mag(i) = abs(output(idx) / input(idx));
end

subplot(223)
plot(freq,phase_shift*180/pi,'.')
ylabel('Phase shift (deg)');
xlabel('Frequency (Hz)');
subplot(224)
plot(freq,mag,'.')
xlabel('Frequency (Hz)');
ylabel('Magnitude');
keyboard

function params = parseVarArgs(params,varargin)
% Parse variable input arguments supplied in name/value format.
%
%    params = parseVarArgs(params,'property1',value1,'property2',value2) sets
%    the fields propertyX in p to valueX.
%
%    params = parseVarArgs(params,varargin{:},'strict') only sets the field
%    names already present in params. All others are ignored.
%
%    params = parseVarArgs(params,varargin{:},'assert') asserts that only
%    parameters are passed that are set in params. Otherwise an error is
%    thrown.
%
% AE 2009-02-24

% check if correct number of inputs
if mod(length(varargin),2)
    switch varargin{end}
        case 'strict'
            % in 'strict' case, remove all fields that are not already in params
            fields = fieldnames(params);
            ndx = find(~ismember(varargin(1:2:end-1),fields));
            varargin([2*ndx-1 2*ndx end]) = [];
        case 'assert'
            % if any fields passed that aren't set in the input structure,
            % throw an error
            extra = setdiff(varargin(1:2:end-1),fieldnames(params));
            assert(isempty(extra),'Invalid parameter: %s!',extra{:})
            varargin(end) = [];
        otherwise
            err.message = 'Name and value input arguments must come in pairs.';
            err.identifier = 'parseVarArgs:wrongInputFormat';
            error(err)
    end
end

% parse arguments
for i = 1:2:length(varargin)
    if ischar(varargin{i})
        params.(varargin{i}) = varargin{i+1};
    else
        err.message = 'Name and value input arguments must come in pairs.';
        err.identifier = 'parseVarArgs:wrongInputFormat';
        error(err)
    end
end