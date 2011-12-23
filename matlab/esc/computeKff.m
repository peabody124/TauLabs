function [esc t rpm setpoint] = computeKff(varargin)
% Computes the feedforward model for your motor

params.esc = [];
params.port = '/dev/tty.usbmodemfd131';
params.dt_ms = 1;
params.low_speed_rpm = 1000;
params.high_speed_rpm = 6000;
params.test_time_s = 15;
params = parseVarArgs(params,varargin{:});

% If no esc object passed in, create one
esc = params.esc;
if isempty(params.esc)
    esc = EscSerial;
end

% Open port if not open
if ~isOpen(esc)
    esc = openPort(esc,params.port);
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
esc = setSerialSpeed(esc, 450);
pause(4);
for i = 450:5:params.low_speed_rpm
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
duty_cycle = [];
low = 0;

tic;
a = [];
while toc < params.test_time_s
    pause(params.dt_ms / 1000);
    speed = params.low_speed_rpm + toc * (params.high_speed_rpm - params.low_speed_rpm) / params.test_time_s;
    esc = setSerialSpeed(esc, speed);
    
    [esc t_ rpm_ setpoint_ dc_] = parseLogging(esc);    
    t = [t t_];
    rpm = [rpm rpm_];
    duty_cycle = [duty_cycle dc_];
    
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

b = robustfit(rpm', duty_cycle')
Kff = b(2) * 32176 / 1024 * 32
Kff2 = -b(1) * 32176 / 1024

esc.configuration.Kff = Kff;
esc.configuration.Kff2 = Kff2;
esc.configuration.RpmMax = (esc.configuration.MaxDc - b(1)) / b(2);
a = [min(rpm) max(rpm)];
plot(rpm,duty_cycle,'.',a,a * b(2) + b(1))

return;
% Run the code below here to test the open loop response
esc.configuration.RisingKp = 0;
esc.configuration.FallingKp = 0;
esc.configuration.Ki = 0;
esc = setConfiguration(esc);

% flush old data and enable logging 
if get(esc.ser,'BytesAvailable') > 0
    fread(esc.ser, get(esc.ser,'BytesAvailable'), 'uint8');
end
esc = enableLogging(esc);
esc.packet = [];

t = [];
rpm = [];
setpoint = [];
duty_cycle = [];
low = 0;

tic;
a = [];
while toc < params.test_time_s
    pause(params.dt_ms / 1000);
    speed = params.low_speed_rpm + toc * (params.high_speed_rpm - params.low_speed_rpm) / params.test_time_s;
    esc = setSerialSpeed(esc, speed);
    
    [esc t_ rpm_ setpoint_ dc_] = parseLogging(esc);    
    t = [t t_];
    rpm = [rpm rpm_];
    duty_cycle = [duty_cycle dc_];
    
    if mean(rpm) < 100
        disp('Failed to start or shut down');
        esc = setSerialSpeed(esc, 0);
        break;
    end
    setpoint = [setpoint setpoint_];
end

keyboard
clf
subplot(211);
plot((t-t(1))/1000,setpoint,(t-t(1))/1000,rpm + 2* params.sin_amp_rpm);
legend('setpoint [rpm]','speed [rpm]')
title('Input and Output (intentionally offset)')
ylim([params.base_speed_rpm-params.sin_amp_rpm*2  params.base_speed_rpm + 4* params.sin_amp_rpm]);
xlabel('Time (s)')
ylabel('RPM');


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