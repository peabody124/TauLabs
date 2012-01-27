
speeds = 600:50:5800

for i = 1:length(speeds)
    esc = setSerialSpeed(esc, speeds(i));
    pause(0.05);
    status = getStatus(esc);
    
    speed(i) = status.CurrentSpeed;
    current(i) = status.Current;
    zcd_fraction(i) = status.ZcdFraction;
end

esc = setSerialSpeed(esc, 600);