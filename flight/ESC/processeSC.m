%dump binary memory dump.bin &back_buf[0] &back_buf[8095]


fid = fopen('dump.bin');
dat = fread(fid);
fclose(fid)
%dat = reshape(dat,2,[]);
dat = double(typecast(uint8(dat(:)),'uint16'));

ref = mean(dat,2);
h(1) = subplot(311)
plot(dat(:,1) - ref,'.');
h(2) = subplot(312);
plot(dat(:,2) - ref,'.');
h(3) = subplot(313);
plot(dat(:,3) - ref,'.');
linkaxes(h,'x');

dat = reshape(dat(1:8092),4,[])';

sync = find(dat == 65535);

state = dat(sync + 1);
detected = dat(sync + 2);
commutated = dat(sync + 3);
dat(sync+1) = 65535;
dat(sync+2) = 65535;
dat(sync+3) = 65535;
dat(sync(end):end) = 65535;

bad = find(mod(diff(sync)-4,3) ~= 0);
for i = 1:length(bad)
    dat(sync(bad(i)):sync(bad(i)+1)) = 65535;
end

dat(dat == 65535) = [];

dat = reshape(dat,3,[]);
plot(1:size(dat,2),dat','.',sync/3,state*40,'k')

sync = sync - 2*(0:length(sync)-1)';

dev = double(dat) - 580;

t = 1:size(dev,2);
h(1) = subplot(211);
plot(t,bsxfun(@minus,bsxfun(@plus,0:2000:4000,dev'),0),sync / 3, state * 100 - 1000, '.k',sync / 3, mod(detected,2) * 1000 - 2000,'.',sync / 3, commutated * 1000 - 2500,'.')
h(2) = subplot(212);
avg = mean(dat,1);
plot(t,bsxfun(@plus,bsxfun(@minus,dat,avg),[0 1500 3000]'),'.');
linkaxes(h,'x')
