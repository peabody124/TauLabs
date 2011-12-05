function [fall_info rise_info] = StepResponse(varargin)
% Computes the step resposnse of an ESC
% Optional params
%  esc - the already existing esc object
%  port - the port to open if there isn't one
%  step_period_ms - the period between sending updates
%  step_size_rpm - the number of rpm to increment by
%  step_base_rpm - the starting prm to ramp up to
%  repetitions - how many steps to do
%
% JC 2011-12-03

params.esc = [];
params.port = '/dev/tty.usbmodemfa141';
params.step_period_ms = 500;
params.step_size_rpm = 500;
params.base_speed_rpm = 2500;
params.repetitions = 20;
params = parseVarArgs(params, varargin{:});

% If no esc object passed in, create one
esc = params.esc;
if isempty(params.esc)
    esc = EscSerial;
end

base_speed_rpm = params.base_speed_rpm;
step_size_rpm = params.step_size_rpm;
step_period_ms = params.step_period_ms;
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
esc = setSerialSpeed(esc, 450);
pause(4);
for i = 550:5:base_speed_rpm
    esc = setSerialSpeed(esc, i);
    pause(0.01);
end

% flush old data and enable logging 
if(get(esc.ser,'BytesAvailable') > 0)
    fread(esc.ser, get(esc.ser,'BytesAvailable'), 'uint8');
end
esc = enableLogging(esc);
esc.packet = [];

clf
subplot(211);
h(1) = line('XData', (1:1),'YData', sin(1:1), 'Color', 'b');
h(2) = line('XData', (1:1),'YData', sin(1:1), 'Color', 'g');
legend('speed [rpm]', 'setpoint [rpm]')

t = [];
rpm = [];
setpoint = [];

low = 0;

tic;
for rep = 1:(params.repetitions*2)
    while toc < (rep * params.step_period_ms / 1000 / 2)
        [esc t_ rpm_ setpoint_] = parseLogging(esc);
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
        ylim([base_speed_rpm - 200 base_speed_rpm + step_size_rpm + 200]);
        xlim([t(max(1,end-4000+1)) t(end)])
        drawnow;
    end
    new_speed = base_speed_rpm + step_size_rpm * (mod(rep,2)==0);
    esc = setSerialSpeed(esc, new_speed);
    disp(['Setting speed to ' num2str(new_speed) ' after ' num2str(toc) ' seconds']);
end

% Shut down controller and resume PWM control.
esc = setSerialSpeed(esc, 0);
esc = disableSerialControl(esc);
esc = disableLogging(esc);
% If we opened the port close it again
if opened
    esc = closePort(esc);
end

step_length = 100;

rise = find(setpoint(1:end-1-step_length) == base_speed_rpm & setpoint(2:end-step_length) == (base_speed_rpm + step_size_rpm));
fall = find(setpoint(1:end-1-step_length) == (base_speed_rpm + step_size_rpm) & setpoint(2:end-step_length) == base_speed_rpm);

rise_response = mean(rpm(bsxfun(@plus,rise,(1:step_length)')),2);
fall_response = mean(rpm(bsxfun(@plus,fall,(1:step_length)')),2);

rise_info = stepinfo(rise_response-base_speed_rpm,t(1:step_length));
fall_info = stepinfo(fall_response-base_speed_rpm-step_size_rpm,t(1:step_length));

disp(['Rise time was: ' num2str(rise_info.RiseTime) ' ms and overshoot was ' num2str(rise_info.Overshoot / step_size_rpm * 100) '%']);
disp(['Fall time was: ' num2str(fall_info.RiseTime) ' ms and overshoot was ' num2str(fall_info.Overshoot / step_size_rpm * 100) '%']);

subplot(212);
plot(t(1:step_length) - t(1), rise_response, t(1:step_length)- t(1), fall_response);


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
