
classdef EscSerial
    properties
        ser = [];
        logging = false;
        port = '';
        configuration = struct( ...
            'RisingKp',int16(0),...
            'FallingKp',int16(0),...
            'Ki',int16(0), ...
            'Kv',int16(850), ...
            'Kff2',int16(0), ...
            'ILim',int16(500), ...
            'MaxError',uint16(250), ...
            'MaxDcChange',int16(100), ...
            'MaxDc',int16(1024 * 0.95), ...
            'MinDc',int16(0), ...
            'InitialStartupSpeed',int16(30), ...
            'FinalStartupSpeed',int16(400), ...
            'StartupCurrentTarget',int16(400), ...
            'SoftCurrentLimit',uint16(40000), ...
            'HardCurrentLimit',uint16(40000), ...
            'CommutationPhase',int16(10), ...
            'CommutationOffset',int16(0), ...
            'PwmMin',uint16(1050), ...
            'PwmMax',uint16(2000), ...
            'RpmMin',uint16(500), ...
            'RpmMax',uint16(10000), ...
            'PwmFreq',uint16(40000), ...
            'Braking',int8(0), ...
            'Direction',int8(0), ...
            'Mode',int8(0));
            
        packet = [];
    end
    
    properties(Constant)
        SYNC_BYTE = hex2dec('85');
        ESC_COMMAND_SET_CONFIG = 0;
        ESC_COMMAND_GET_CONFIG = 1;
        ESC_COMMAND_SAVE_CONFIG = 2;
        ESC_COMMAND_ENABLE_SERIAL_LOGGING = 3;
        ESC_COMMAND_DISABLE_SERIAL_LOGGING = 4;
        ESC_COMMAND_REBOOT_BL = 5;
        ESC_COMMAND_ENABLE_SERIAL_CONTROL = 6;
        ESC_COMMAND_DISABLE_SERIAL_CONTROL = 7;
        ESC_COMMAND_SET_SPEED = 8;
        ESC_COMMAND_WHOAMI = 9;
        ESC_COMMAND_ENABLE_ADC_LOGGING = 10;
        ESC_COMMAND_GET_ADC_LOG = 11;
        ESC_COMMAND_SET_PWM_FREQ = 12;
        ESC_COMMAND_GET_STATUS = 13;
        ESC_COMMAND_BOOTLOADER = 14;
        ESC_COMMAND_GETVOLTAGES = 15;
        ESC_COMMAND_LAST = 16;
        
        READ_PACKET_LENGTH = 6*2;
    end
    
    methods
        function self = EscSerial()
        end
        
        function self = openPort(self,serport)
            self.ser = serial(serport);
            self.ser.BaudRate=115200; %230400;
            self.ser.InputBufferSize=2*1024*1024;
            self.port = serport;
            fopen(self.ser);
        end
        
        function self = closePort(self)
            fclose(self.ser)
            delete(self.ser)
            self.ser = [];
        end
        
        function open = isOpen(self)            
            open = ~isempty(self.ser) && all(get(self.ser,'status')=='open');
        end
        
        function self = enableLogging(self)
            assert(isOpen(self), 'Open serial port first');
            fwrite(self.ser, uint8([self.SYNC_BYTE self.ESC_COMMAND_ENABLE_SERIAL_LOGGING]));
            % TODO: Get/check ack
            self.logging = true;
            self.packet = [];
        end
        
        function self = disableLogging(self)
            assert(isOpen(self), 'Open serial port first');
            self.logging = false;
            
            fwrite(self.ser, uint8([self.SYNC_BYTE self.ESC_COMMAND_DISABLE_SERIAL_LOGGING]));
            % TODO: Get/check ack
            pause(0.1);
            %fread(self.ser,'uint8'); % throw away remaining data
        end
        
        function self = enableSerialControl(self)
            assert(isOpen(self), 'Open serial port first');
            fwrite(self.ser, uint8([self.SYNC_BYTE self.ESC_COMMAND_ENABLE_SERIAL_CONTROL]));
            % TODO: Get/check ack
        end
        
        function self = disableSerialControl(self)
            assert(isOpen(self), 'Open serial port first');
            fwrite(self.ser, uint8([self.SYNC_BYTE self.ESC_COMMAND_DISABLE_SERIAL_CONTROL]));
            % TODO: Get/check ack
        end

        function self = setSerialSpeed(self, speed)
            assert(isOpen(self), 'Open serial port first');
            fwrite(self.ser, uint8([self.SYNC_BYTE self.ESC_COMMAND_SET_SPEED typecast(int16(speed), 'uint8')]));
            % TODO: Get/check ack
        end
        
        function self = setConfiguration(self)
            assert(isOpen(self), 'Open serial port first');
            assert(self.configuration.MaxDc <= 1024, 'MaxDc is over the safe limit');
            command = [uint8([self.SYNC_BYTE self.ESC_COMMAND_SET_CONFIG]), ...
                typecast(int16(self.configuration.RisingKp),'uint8'), ...
                typecast(int16(self.configuration.FallingKp),'uint8'), ...
                typecast(int16(self.configuration.Ki),'uint8'), ...
                typecast(int16(self.configuration.Kv),'uint8'), ...
                typecast(int16(self.configuration.Kff2),'uint8'), ...
                typecast(int16(self.configuration.ILim),'uint8'), ...
                typecast(uint16(self.configuration.MaxError),'uint8'), ...  
                typecast(int16(self.configuration.MaxDcChange),'uint8'), ...
                typecast(int16(self.configuration.MaxDc),'uint8'), ...
                typecast(int16(self.configuration.MinDc),'uint8'), ...
                typecast(uint16(self.configuration.InitialStartupSpeed),'uint8'), ...
                typecast(uint16(self.configuration.FinalStartupSpeed),'uint8'), ...
                typecast(uint16(self.configuration.StartupCurrentTarget),'uint8'), ...
                typecast(uint16(self.configuration.SoftCurrentLimit),'uint8'), ...
                typecast(uint16(self.configuration.HardCurrentLimit),'uint8'), ...
                typecast(int16(self.configuration.CommutationPhase),'uint8'), ...
                typecast(int16(self.configuration.CommutationOffset),'uint8'), ...
                typecast(uint16(self.configuration.PwmMin),'uint8'), ...
                typecast(uint16(self.configuration.PwmMax),'uint8'), ...
                typecast(uint16(self.configuration.RpmMin),'uint8'), ...
                typecast(uint16(self.configuration.RpmMax),'uint8'), ...
                typecast(uint16(self.configuration.PwmFreq),'uint8'), ...
                typecast(uint8(self.configuration.Braking),'uint8'), ...
                typecast(uint8(self.configuration.Direction),'uint8'), ...
                typecast(uint8(self.configuration.Mode),'uint8'), ...
                ];
            % TODO: Receive/Check ack
            fwrite(self.ser, command(1:15));
            pause(0.1)
            fwrite(self.ser, command(16:end));
        end
        
        function config = getConfiguration(self)
            % Get the configuration from the ESC
            % dat = getConfiguration(self)
            assert(isOpen(self), 'Open serial port first');
            assert(~self.logging, 'Do not do this while serial logging is enabled');
            flush(self);
            command = uint8([self.SYNC_BYTE self.ESC_COMMAND_GET_CONFIG]);
            fwrite(self.ser, command);
            pause(0.1);
            dat = uint8(fread(self.ser, 47, 'uint8'));
            config.RisingKp = typecast(dat(1:2),'int16');
            config.FallingKp = typecast(dat(3:4),'int16');
            config.Ki = typecast(dat(5:6),'int16');
            config.Kv = typecast(dat(7:8),'int16');
            config.Kff2 = typecast(dat(9:10),'int16');
            config.ILim = typecast(dat(11:12),'int16');
            config.MaxError = typecast(dat(13:14),'uint16');
            config.MaxDcChange = typecast(dat(15:16),'int16');
            config.MaxDc = typecast(dat(17:18),'int16');
            config.MinDc = typecast(dat(19:20),'int16');
            config.InitialStartupSpeed = typecast(dat(21:22),'uint16');
            config.FinalStartupSpeed = typecast(dat(23:24),'uint16');
            config.StartupCurrentTarget = typecast(dat(25:26),'uint16');
            config.SoftCurrentLimit = typecast(dat(27:28),'uint16');
            config.HardCurrentLimit = typecast(dat(29:30),'uint16');
            config.CommutationPhase = typecast(dat(31:32),'int16');
            config.CommutationOffset = typecast(dat(33:34),'int16');
            config.PwmMin = typecast(dat(35:36),'uint16');
            config.PwmMax = typecast(dat(37:38),'uint16');
            config.RpmMin = typecast(dat(39:40),'uint16');
            config.RpmMax = typecast(dat(41:42),'uint16');
            config.PwmFreq = typecast(dat(43:44),'uint16');
            config.Braking = typecast(dat(45), 'uint8');
            config.Direction = typecast(dat(46),'uint8');
            config.Mode = typecast(dat(47),'uint8');
        end
        
        function self = saveConfiguration(self)
            assert(isOpen(self), 'Open serial port first');
            fwrite(self.ser, uint8([self.SYNC_BYTE self.ESC_COMMAND_SAVE_CONFIG]));
            retval = fread(self.ser, 1, 'int8')
        end
        
        function status = getStatus(self)
            % Get the configuration from the ESC
            % dat = getConfiguration(self)
            assert(isOpen(self), 'Open serial port first');
            assert(~self.logging, 'Do not do this while serial logging is enabled');
            flush(self);
            command = uint8([self.SYNC_BYTE self.ESC_COMMAND_GET_STATUS]);
            fwrite(self.ser, command);
            pause(0.1);
            dat = uint8(fread(self.ser, 17, 'uint8'));
            status.SpeedSetpoint = typecast(dat(1:2),'int16');
            status.CurrentSpeed = typecast(dat(3:4),'uint16');
            status.Current = typecast(dat(5:6),'uint16');
            status.TotalCurrent = typecast(dat(7:8),'uint16');
            status.Kv = typecast(dat(9:10),'int16');
            status.Battery = typecast(dat(11:12),'uint16');
            status.TotalMissed = typecast(dat(13:14),'uint16');
            status.DutyCycle = dat(15);
            status.ZcdFraction = dat(16);
            status.Error = dat(17);
        end
        
        function voltages = getVoltages(self)
            % Get the configuration from the ESC
            % dat = getConfiguration(self)
            assert(isOpen(self), 'Open serial port first');
            assert(~self.logging, 'Do not do this while serial logging is enabled');
            flush(self);
            command = uint8([self.SYNC_BYTE self.ESC_COMMAND_GETVOLTAGES]);
            fwrite(self.ser, command);
            pause(0.1);
            voltages = fread(self.ser, 3*6, 'int16');
            voltages = reshape(voltages,3,6);
        end
        
        function [self t rpm setpoint dutycycle] = parseLogging(self)
            assert(isOpen(self), 'Open serial port first');
            assert(self.logging, 'Enable logging first');
            
            count = get(self.ser, 'BytesAvailable');
            if  count> 0
                self.packet = [self.packet; uint8(fread(self.ser, count,'uint8'))];
            end
            
            head = find(self.packet(1:end-1) == 0 & self.packet(2:end) == 255, 1, 'first');
            j = 1;
            t = [];
            rpm = [];
            setpoint = [];
            dutycycle = [];
            
            while(~isempty(head) && head < (length(self.packet) - self.READ_PACKET_LENGTH))
                
           		a = uint8((self.packet(head+2:head+5)));
                t(j) = double(typecast(a,'uint32'));
				
                a = uint8((self.packet(head+6:head+7)));
                rpm(j) = double(typecast(a,'uint16'));
		
                a = uint8((self.packet(head+8:head+9)));
                setpoint(j) = double(typecast(a,'int16'));
		
                a = uint8((self.packet(head+10:head+11)));
                dutycycle(j) = double(typecast(a,'int16'));
                

                self.packet(1:head+11) = [];
                head = find(self.packet(1:end-1) == 0 & self.packet(2:end) == 255, 1, 'first');
                j = j + 1;
            end 
        end
        
        function self = flush(self)
            % flush old data
            if get(self.ser, 'BytesAvailable') > 0
                fread(self.ser, get(self.ser,'BytesAvailable'), 'uint8');
            end
        end
        
        function self = plotLogging(self)

            assert(get(self.ser,'BytesAvailable') > 0, 'Logging may not be enabled as no bytes available');
            
            % flush old data
            fread(self.ser, get(self.ser,'BytesAvailable'), 'uint8');

            h(1) = line('XData', (1:1),'YData', sin(1:1), 'Color', 'b');
            h(2) = line('XData', (1:1),'YData', sin(1:1), 'Color', 'g');
            legend('speed [rpm]', 'setpoint [rpm]')

            t = [];
            rpm = [];
            setpoint = [];
            self.packet = [];
            for i = 1:100
                [self t_ rpm_ setpoint_] = parseLogging(self);
                t = [t t_];
                rpm = [rpm rpm_];
                setpoint = [setpoint setpoint_];

                
                set(h(1), 'XData', t(max(1,end-4000):end), 'YData', rpm(max(1,end-4000):end))
                set(h(2), 'XData', t(max(1,end-4000):end), 'YData', setpoint(max(1,end-4000):end))
                xlim([t(max(1,end-4000+1)) t(end)])
                drawnow;

                
                shg
            end
        end
        
        function self = setPwmRate(self, freq)
            assert(isOpen(self), 'Open serial port first');
            fwrite(self.ser, uint8([self.SYNC_BYTE self.ESC_COMMAND_SET_PWM_FREQ typecast(uint16(freq), 'uint8')]));
            % TODO: Get/check ack
        end
        
        function [self dat] = runAdcLog(self)
            assert(isOpen(self), 'Open serial port first');
            assert(~self.logging, 'Do not run this while doing serial logging');
            fwrite(self.ser, uint8([self.SYNC_BYTE self.ESC_COMMAND_ENABLE_ADC_LOGGING]));
            flush(self);
            % TODO: Get/check ack
            pause(0.5);
            fwrite(self.ser, uint8([self.SYNC_BYTE self.ESC_COMMAND_GET_ADC_LOG]));
            for i = 1:40
                if get(self.ser,'BytesAvailable') == 12048
                    break;
                end
                pause(0.1);
            end
            assert(i~=100,'Timed out getting ADC data');
            dat = uint8(fread(self.ser, get(self.ser,'BytesAvailable'), 'uint8'));
            dat = typecast(dat,'uint16');
            dat = double(reshape(dat,4,[]));
        end
        
        function self = reprogram(self)
            assert(isOpen(self), 'Open serial port first');
            fwrite(self.ser, uint8([self.SYNC_BYTE self.ESC_COMMAND_BOOTLOADER hex2dec('73') hex2dec('37')]));
            self = closePort(self);
            escbin;
            system(['./tools/stm32flash/stm32flash -w esc.bin -v -g 0x0 -b 115200 ' self.port]);
            self = openPort(self, self.port);
        end

        function self = program(self)
            self = closePort(self);
            escbin;
            system(['./tools/stm32flash/stm32flash -w esc.bin -v -g 0x0 -b 115200 ' self.port]);
            self = openPort(self, self.port);
        end

        function self = boot(self)
            self = closePort(self);
            escbin;
            system(['./tools/stm32flash/stm32flash -g 0x0 -b 115200 ' self.port]);
            self = openPort(self, self.port);
        end

        function display(esc)
            if(isOpen(esc))
                disp('Esc is connected');
            else
                disp('Esc is not connected');
            end
        end
    end
end

            
