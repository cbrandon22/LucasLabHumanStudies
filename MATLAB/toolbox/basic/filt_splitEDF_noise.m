function filt_splitEDF_noise
fn = '/Volumes/HumanStudies/HumanStudies/CCDT/eeg/HUP168/eeg.noreref/HUP168_23Apr18_1051.258';
srate = 1024;

c = strsplit(fn,'.');
splitfn = c(1:end-1);
fnbase = splitfn{1};
for i=2:length(splitfn)
    fnbase=[fnbase '.' splitfn{i}];
end
c = str2num(c{end});
dat = look(fnbase,c,0);
ln = mtmlinenoise(dat,3,srate,srate,60:60:300);% get line noise
dat2 = dat-ln;
% figure;
% plot(dat2)
dat3=int16(dat2);
fchan = fopen(fn,'w','l');
fwrite(fchan,dat3,'int16');
fclose(fchan);
%dat4 = look(fnbase,258);