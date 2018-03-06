function data = loadtxtfile(fname)
% function data = loadtxtfile(fname)
%   Load data saved to text file by wireless recorder.

% parameters
voltage_reference = 1.2;
adc_resolution = 8;
input_prescaling = 1/3;
AFE_gain = 2000;

% load data
fid = fopen(fname,'r');
data = [];
for i = 1:4 % tank header
    fgetl(fid);
end
while ~feof(fid)
    for i = 1:3 % block header
        fgetl(fid);
    end
    bk = fscanf(fid,'%x',[16 16]); % read block (16 bit)
    if all(size(bk)==[16 16]), data = [data, bk]; end
end
fclose(fid);
data = reshape(data,1,size(data,2)*16);
data = (data/2^adc_resolution)*(voltage_reference/AFE_gain)*(1/input_prescaling)*1e6; % microvolts
