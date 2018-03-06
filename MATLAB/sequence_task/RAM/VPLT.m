function VPLT
% function VPLT
%   Visual stimulus locations and sizes assume a viewing distance of ~ 30
%   inches to get the correct visual angle: VA = 2*atan(size/(2*30)).
%
%   DR 12/2014

% parameters
stimset = 50; % stim set number
stimdir = 'C:\TDT\Matlab\VPLT\stimsets\'; % stim set directory
vt_size = [25 25]; % visual target size (for CCDT)
vt_cfix = ntsc2rgb([0.8 0 0]); % visual target fixation color
vt_cres = ntsc2rgb([0.8 .3 -.3]); % visual target response color
pd_cfix = [0.7 0.7 0.7]; % photodetector target fixation color
pd_cres = [1.0 1.0 1.0]; % photodetector target response color
usb_bardr = 19200; % serial port baud rate
usb_ptnum = 7; % serial port number (COM-X in Device Manager)

% load all 200 images into memory
cd([stimdir 'SET' sprintf('%03d',stimset)]); % open stim set directory
C = imread('1.bmp'); Cblack = zeros(size(C),'uint8'); Cpict = zeros([size(C),200],'uint8'); % assumes 200 images
for ii = 1:200
    C = imread([num2str(ii) '.bmp']);
    Cpict(:,:,:,ii) = reshape(C,[size(C),1]);
end

% create figure
hf = figure('Name','','NumberTitle','off','Toolbar','none','Menubar','none',...
            'Resize','off','DoubleBuffer','off','Renderer','painters','Units','pixel',...
            'Position',[-1919 -29 1920 1080],'Color','k');
himg = image(Cblack); axis square;
set(gca,'Units','pixels','Position',[533 281 923 535],'Visible','off','XLimMode','manual',... % ~11 deg VA
        'YLimMode','manual','ZLimMode','manual','CLimMode','manual','ALimMode','manual');
hvt(1) = axes('Parent',hf,'Units','pixel','Position',[727, 281, vt_size],'Color',vt_cfix,'Visible','off'); % ~0.5 deg VA
hvt(2) = axes('Parent',hf,'Units','pixel','Position',[727+535-vt_size(1)-1, 281, vt_size],'Color',vt_cfix,'Visible','off');
hvt(3) = axes('Parent',hf,'Units','pixel','Position',[727, 281+535-vt_size(2)-1, vt_size],'Color',vt_cfix,'Visible','off');
hvt(4) = axes('Parent',hf,'Units','pixel','Position',[727+535-vt_size(1)-1, 281+535-vt_size(2)-1, vt_size],'Color',vt_cfix,'Visible','off');
hvt(5) = axes('Parent',hf,'Units','pixel','Position',[727+535/2-vt_size(1)/2, 281+535/2-vt_size(2)/2, vt_size],'Color',vt_cfix,'Visible','off');
hfx = axes('Parent',hf,'Units','pixel','Position',[727+535/2-25, 281+535/2-25, 50, 50],'Color',vt_cfix,'Visible','off'); % ~1 deg VA
hpd = axes('Parent',hf,'Units','pixel','Position',[0 0 100 90],'Color',pd_cfix,'Visible','off');
set(hpd,'XLimMode','manual','YLimMode','manual','ZLimMode','manual',...
        'CLimMode','manual','ALimMode','manual');
drawnow;

% open serial port
usbpt = serial(['com' num2str(usb_ptnum)],'baudrate',usb_bardr,'parity','none','databits',8);
usbpt.ReadAsyncMode = 'continuous';
fopen(usbpt);

% connect to port and interpret 8-bit code
% 0 = VPLT intertrial
% 1-200 = VPLT display image 1-200
% 201 = VPLT fixation
% 210 = CCDT intertrial
% 211-215 = CCDT fixation
% 216 = CCDT response
% 217 or greater = stop task
lastCMD = 0;
fstop = 0;
while ~fstop
    if usbpt.BytesAvailable
        out = fread(usbpt,usbpt.BytesAvailable,'char'); % 8-bit character
        nowCMD = out(end);
        if nowCMD~=lastCMD
            if (nowCMD>0) && (nowCMD<=200) % VPLT: display image
                set(hfx,'Visible','off');
                set(himg,'CData',squeeze(Cpict(:,:,:,nowCMD)));
                set(hpd,'Color',pd_cres);
            elseif nowCMD==0 % VPLT: intertrial
                set(himg,'CData',Cblack);
                set(hpd,'Visible','off');
                set(hpd,'Color',pd_cfix);
            elseif nowCMD==201 % VPLT: fixation
                set(hfx,'Visible','on');
                set(hpd,'Visible','on');
            elseif nowCMD==210 % CCDT: intertrial
                set(hvt,'Visible','off');
                set(hpd,'Visible','off');
                set(hvt,'Color',vt_cfix);
                set(hpd,'Color',pd_cfix);
            elseif (nowCMD>210) && (nowCMD<=215) % CCDT: fixation
                set(hvt(nowCMD-210),'Visible','on');
                set(hpd,'Visible','on');
            elseif nowCMD==216 % CCDT: response
                set(hvt,'Color',vt_cres);
                set(hpd,'Color',pd_cres);
            elseif nowCMD==217
                fstop = 1;
            end
            drawnow
            lastCMD = nowCMD;
        end
    end
end

% close serial port and figure
fclose(usbpt);
close(gcf);