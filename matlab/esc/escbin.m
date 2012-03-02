bl_size = hex2dec('2000');
fid = fopen('build/bl_esc/bl_esc.bin');
dat = uint8(fread(fid,'uchar'));
fclose(fid);
dat(bl_size) = 0;  % pad till end
fid = fopen('build/fw_esc/fw_esc.bin');
dat = [dat; uint8(fread(fid,'uchar'))];
fclose(fid);
fid = fopen('esc.bin','w');
fwrite(fid,dat)
fclose(fid)