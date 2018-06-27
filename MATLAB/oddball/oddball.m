function oddball
%% LabJack setup
% Creates labjack object using labjack.m. Set silentMode to true in object
% construction ('silentMode',true) to run without communicating with LabJack or saving TTL/eye
% tracker data. Automatically enters silentMode if no LabJack connected.
events = {'BACKGROUNDLF','BACKGROUNDHF','TARGETLF','TARGETHF','GO','NO_GO',...
    'RESPONSE','PASSIVE_START','PASSIVE_END','RESPONSE_START','RESPONSE_END',...
    'START_INSTRUCTIONS','DONE_INSTRUCTIONS','SESS_START','SESS_END'}
daq = labJack('verbose',true);

% LabJack FIO port settings
daqinfo.pulsePort = 0; %TTL pulse port output
daqinfo.xPort = 4; %x coordinate voltage input
daqinfo.yPort = 5; %y coordinate voltage input
daqinfo.diamPort = 6; %pupil diameter voltage input
daqinfo.pulsewidth = 5; %milliseconds
daqinfo.recBuff = zeros(16,1); % Read bytes (length = 10+ 2*Num of channels to read AIN)

% ConfigIO byte. See Labjack U3 documentation (5.2.3)
% Sets all channels other than FIO0 and FIO1 to analog input, Use ioBuffer[11] to change channel
% config (8 bit command, 0=digital,1=analog ex: 252 = 11111100)
ioBuffer = [0 248 3 11 0 0 13 0 64 0 252 0]';
ioBuffer = daq.checksum(ioBuffer,'extended');
bytesWritten = daq.rawWrite(ioBuffer);
if bytesWritten == 0, warning('ConfigIO error: No bytes written.');end
if bytesWritten < 12, warning('ConfigIO error: Did not write all of the buffer.');end
recBuffer = zeros(12,1);
[bytesRead,retBytes] = daq.rawRead(recBuffer,12);
if bytesRead == 0, warning('ConfigIO error: Read failed.');end
if bytesRead < 12, warning('ConfigIO error: Did not read all of the buffer.');end

% Bytes to send to read AIN from eye tracker channels
sendBuff = [0 248 5 0 0 0 0 1 daqinfo.xPort 31 1 daqinfo.yPort 31 1 daqinfo.diamPort 31]';
sendBuff = daq.checksum(sendBuff,'extended');


function sendPulse
    % Send TTL pulse
    if daq.silentMode, return; end
    time = GetSecs; %getting timestamp. Cmd-response for labjack is <1ms
    daq.timedTTL(daqinfo.pulsePort, daqinfo.pulsewidth);
    syncTimes = [syncTimes; time];


function readInput(target)
    % Writes to then reads from labjack for eyetracker input
    if daq.silentMode, return; end
    bytesWritten = daq.rawWrite(sendBuff);
    if bytesWritten == 0, warning('Feedback error: No bytes written.');end
    if bytesWritten < 16, warning('Feedback error: Did not write all of the buffer.');end
    [bytesRead, retBytes] = daq.rawRead(daqinfo.recBuff,16);
    if bytesRead == 0, warning('Feedback error: Read failed.');end
    if bytesRead < length(daqinfo.recBuff), warning('Feedback error: Did not read all of the buffer.');end
    xVolts = double(retBytes(11))*256 + double(retBytes(10));
    yVolts = double(retBytes(13))*256 + double(retBytes(12));
    diamVolts = double(retBytes(15))*256 + double(retBytes(14));
    eyeData = [eyeData;[xVolts,yVolts,diamVolts,target]];
