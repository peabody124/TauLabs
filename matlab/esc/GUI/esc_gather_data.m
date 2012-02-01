function esc_gather_data(handles)

cla

if (0) %Dummy function which just plots lines
	%Initialize line arrays
	h_hallSensors = line('XData', (1),'YData', sin(1), 'Color', 'b','LineStyle','o');
	h_fluxEstimation = line('XData', (1),'YData', sin(1), 'Color', 'r','LineStyle','--');

	%DON'T SET LEGEND, IT REALLY SLOWS DOWN THE PLOTTING
% 	legend('RPM (Hall Sensors)', 'RPM (Flux Estimator)', 'Location', 'SE')
	
	%%
	
	startTime=clock;
	
	for i=1:100
		while etime(clock,startTime) < 0.01*i
			drawnow
		end
		set(h_hallSensors, 'XData', (1:0.01:10)+i/10, 'YData', cos((1:0.01:10)+i/10));
		set(h_fluxEstimation, 'XData', (1:0.01:10)+i/10, 'YData', sin((1:0.01:10)+i/10));
		
		set(handles.textFluxEst, 'String', num2str(cos(10+i), 4)); %Update flux estimation speed
		set(handles.textHallSensors, 'String', num2str(sin(10+i), 4)); %Update hall speed
		set(handles.textDataTime, 'String', num2str(10+i, '%0.2f')); %Update time

	end
	
	etime(clock,startTime)
	
else % Function which reads in data from serial port and then plots to 
	  % screen. It's intended for OPESC protocol, so not much use for MIT 
	  % Scooter. However, it would be very useful if the parseLogging
	  % function were made to read and return the serial port data, in a
	  % format that fits into the general structure outlines below

  	global esc

	%Initialize line arrays
	h_speed = line('XData', (1),'YData', sin(1), 'Color', 'b','LineStyle','o');
	h_setpoint = line('XData', (1),'YData', sin(1), 'Color', 'r','LineStyle','--');

	
	if get(handles.checkboxEnableLogging, 'Value')==1
		logOnOff=1;
	else
		logOnOff=0;
	end

	%Allocate storage arrays
	t=0;
	rpm=0;
	duty_cycle=0;
	setpoint=0;
	
	%Initialize storage array indices
	tIdx=1;
	rpmIdx=1;
	duty_cycleIdx=1;
	setpointIdx=1;
	
	%Run measurement loop
	startTime=clock;
	loopCnt=1;
	tic
% 	while(etime(clock,startTime) < 2)
	while(get(handles.buttonRun, 'Value'))
		if toc > .1
			tic
			drawnow
		end
		
% 		while etime(clock,startTime) < 0.5*loopCnt
% 			drawnow
% 			while etime(clock,startTime) < 0.5*loopCnt
% 			end
% 		end
% % 		display([ num2str(etime(clock,startTime)) ' ' num2str(toc) ' ' num2str(0.1*loopCnt)])
		
		% Read serial port logging data
		[esc t_ rpm_ setpoint_ duty_cycle_] = parseLogging(esc);
		if isempty(t_)
			continue;
		end
		
		%------------------------------------------------%
		% Increase array memory allocation, if necessary %
		%------------------------------------------------%
		while length(t) < tIdx+length(t_)
			t(end*2)=0;
		end
		
		while length(rpm) < rpmIdx+length(rpm_)
			rpm(end*2)=0;
		end
		
		while length(duty_cycle) < duty_cycleIdx+length(duty_cycle_)
			duty_cycle(end*2)=0;
		end
		
		while length(setpoint) < setpointIdx+length(setpoint_)
			setpoint(end*2)=0;
		end
		
		%-----------------------------%
		% Assign new values to arrays %
		%-----------------------------%
		t(tIdx:tIdx+length(t_)-1) = t_;
		tIdx=tIdx+length(t_);
		
		rpm(rpmIdx:rpmIdx+length(rpm_)-1) = rpm_;
		rpmIdx=rpmIdx+length(rpm_);
		
		duty_cycle(duty_cycleIdx:duty_cycleIdx+length(duty_cycle_)-1) = duty_cycle_;
		duty_cycleIdx=duty_cycleIdx+length(duty_cycle_);
		
		setpoint(setpointIdx:setpointIdx+length(setpoint_)-1) = setpoint_;
		setpointIdx=setpointIdx+length(setpoint_);
		
		time_slice=max(1,tIdx-5000):2:tIdx-1;
		time_slice=time_slice(time_slice<=time_slice(end));
		%Update GUI
		set(h_speed, 'XData', t(time_slice), 'YData', rpm(time_slice));
		set(h_setpoint, 'XData', t(time_slice), 'YData', setpoint(time_slice));
		
		set(handles.textHallSensors, 'String', num2str(mean(rpm_(max(1,end-30):end)), 4)); %Update hall speed
		set(handles.textFluxEst, 'String', num2str(mean(setpoint_(max(1,end-30):end)))); %Update flux estimation speed
		set(handles.textDataTime, 'String', num2str(t_(end)/1000, '%0.2f')); %Update time

		loopCnt=loopCnt+1;
	end
	
	t=t(1:tIdx-1);
	rpm=rpm(1:rpmIdx-1);
	setpoint=setpoint(1:setpointIdx-1);
	set(h_speed, 'XData', t, 'YData', rpm);
	set(h_setpoint, 'XData', t, 'YData', setpoint);
	
	if logOnOff %If logging is turned on, save the file to a *.mat
		%Get the path destination
		pathstr=fileparts(get(handles.editLogFileName, 'String'));
		
		%Generate log name
		time=clock;
		logName=['escLog ' num2str(time(1), '%04d') '-' num2str(time(2), '%02d') '-' num2str(time(3), '%02d') '-' num2str(time(4), '%02d') ':' num2str(time(5), '%02d') ':' num2str(round(time(6)), '%02d')];
		
		%Save logfile
		save(fullfile(pathstr, logName), 't', 'rpm', 'duty_cycle', 'setpoint');
	end
end