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
esc = setSerialSpeed(esc,0)
esc = disableSerialControl(esc)