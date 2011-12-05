
classdef EscSerial
    properties
        ser = [];
        logging = false;
        configuration = struct( ...
            'RisingKp',int16(40),...
            'FallingKp',int16(100),...
            'Ki',int16(1), ...
            'Kff',int16(3), ...
            'Kff2',int16(0), ...
            'ILim',int16(500), ...
            'MaxDcChange',int16(20), ...
            'MaxDc',int16(321768), ...
            'MinDc',int16(0), ...
            'InitialStartupSpeed',int16(100), ...
            'FinalStartupSpeed',int16(400), ...
            'StartupCurrentTarget',int16(400), ...
            'SoftCurrentLimit',int16(4000), ...
            'HardCurrentLimit',int16(4000), ...
            'CommutationPhase',int16(23), ...
            'CommutationOffset',int16(0), ...
            'Direction',int16(0));
            
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
        ESC_COMMAND_LAST = 13;
        
        READ_PACKET_LENGTH = 6*2;
    end
    
    methods
        function self = EscSerial()
        end
        
        function self = openPort(self,serport)
            self.ser = serial(serport);
            self.ser.BaudRate=115200; %230400;
            self.ser.InputBufferSize=2*1024*1024;
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
            command = [uint8([self.SYNC_BYTE self.ESC_COMMAND_SET_CONFIG]), ...
                typecast(int16(self.configuration.RisingKp),'uint8'), ...
                typecast(int16(self.configuration.FallingKp),'uint8'), ...
                typecast(int16(self.configuration.Ki),'uint8'), ...
                typecast(int16(self.configuration.Kff),'uint8'), ...
                typecast(int16(self.configuration.Kff2),'uint8'), ...
                typecast(int16(self.configuration.ILim),'uint8'), ...
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
                typecast(uint8(self.configuration.Direction),'uint8'), ...
                ];
            % TODO: Receive/Check ack
            fwrite(self.ser, command(1:15));
            pause(0.1)
            fwrite(self.ser, command(16:end));
        end
        
        function dat = getConfiguration(self)
            % Get the configuration from the ESC
            % dat = getConfiguration(self)
            assert(isOpen(self), 'Open serial port first');
            assert(~self.logging, 'Do not do this while serial logging is enabled');
            flush(self);
            command = uint8([self.SYNC_BYTE self.ESC_COMMAND_SET_CONFIG]);
            fwrite(self.ser, command);
            pause(0.1);
            dat = uint8(fread(self.ser, 33, 'uint8'));
        end
        
        function self = saveConfiguration(self)
            assert(isOpen(self), 'Open serial port first');
            fwrite(self.ser, uint8([self.SYNC_BYTE self.ESC_COMMAND_SAVE_CONFIG]));
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
            fwrite(self.ser, uint8([self.SYNC_BYTE self.ESC_COMMAND_SET_PWM_FREQ typecast(int16(freq), 'uint8')]));
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
                if get(self.ser,'BytesAvailable') == 4044
                    break;
                end
                pause(0.1);
            end
            assert(i~=100,'Timed out getting ADC data');
            dat = uint8(fread(self.ser, get(self.ser,'BytesAvailable'), 'uint8'));
            dat = typecast(dat,'uint16');
            dat = double(reshape(dat,3,[]));
        end

    end
end

            
