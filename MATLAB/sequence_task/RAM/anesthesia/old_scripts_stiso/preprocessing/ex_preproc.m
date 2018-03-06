%% Extracting Irvine data
file_path = ['/Volumes/Keep_It_Stored_Keep_It_Secert/Berkeley lab/MNS hands ECoG/IR47_withBallandPassive/raw/2016092015_0001.edf'];
channels_to_save = [1:81];
[header, data] = extract_nk_converted_edf(file_path, channels_to_save);
%[header, data] = extract_nk_converted_edf(file_path, channels_to_save, resample_rate)

%% for Anat's computer:
DTfolderAll = ['/Volumes/Keep_It_Stored_Keep_It_Secert/Berkeley lab/MNS hands ECoG/IR47_withBallandPassive/Passive'];

subj = 'IR47'; block = 'firstRight';
%blocks = [{'raw'}];
Enum = 1:81;

noisy_channels = [70:81]; %epileptic / just don't include in CAR but we will plot
empty_channels = [13:28 69]; %either loose or not recording, I added here HH and AMG channels as wanted to include on grid elecs in CAR 
cut_channels =[1:4]; % These are DC channels we don't want in CAR
reference = [69]; % 1 in Irvine
bad_elecs = unique([noisy_channels empty_channels cut_channels reference]);
good_elecs = setxor(Enum,bad_elecs);

save ([DTfolderAll '/good_elecs'], 'good_elecs')
%to be excluded from CAR
excludeCAR = []; 

%% Create vars
srate = 5000; %in Irvine
ANsrate = 5000;



% create_subj_globals(subj, block,srate, ANsrate, good, bad_elecs, DTfolder_all);
% add_subj_globals(subj, block,'blocks',blocks);
% get_subj_globals('SFC1','firstRight')

%% Create electrode matrix gdat

gdat = data(1:81, :);
ANdata = data(1:4, :);

save([DTfolderAll,'/gdat'],'gdat','-v7.3');
save([DTfolderAll,'/ANdata'],'ANdata','-v7.3');

%% resample
[p, q] = rat(1000/srate); %rational fraction for resampling
%gdat = zeros(numel(Enum),ceil(size(temp_data,2)/(q/p)),'single');
eval('gdat = resample(double(transpose(gdat)), p, q);');
gdat= transpose(gdat);
save([DTfolderAll,'/gdat'],'gdat','-v7.3');

srate = 1000;
save ([DTfolderAll,'/srate'],'srate');
% add_subj_globals('GP44', 'decision', 'srate', srate) % AA- Do we have anymore globals to add?

%% Sanity check of data and srate with FFT plot
%sets up a structure called ecog
ecog.data = gdat;
ecog.srate = srate;

% here I included just the first channel which makes the estimation of the spectrum faster
ecog.data = ecog.data(30,: );


% running the spectral analysis
tic;[S,F] = ecogSpectrum( ecog );toc;
x = squeeze( mean( S,3 ));

% plotting the results
plot( F,x )
%xlim([0,250])

%% Filtering
% Filter LF (<0.1) & HF (>180) noise (noisy data)

[a,b] = butter( 4,[.5/(srate/2), 180/(srate/2)]);
gdat = permute( filtfilt( a,b,permute( double( gdat ),[2 1])),[2 1]);

% Filter 60 Hz & its harmonics

[a,b] = butter( 4,[59/(srate/2), 61/(srate/2)],'stop');
gdat = permute( filtfilt( a,b,permute( double( gdat ),[2 1])),[2 1]);

[a,b] = butter( 4,[119/(srate/2), 121/(srate/2)],'stop');
gdat = permute( filtfilt( a,b,permute( double( gdat ),[2 1])),[2 1]);
 
[a,b] = butter( 4,[179/(srate/2), 181/(srate/2)],'stop');
gdat = permute( filtfilt( a,b,permute( double( gdat ),[2 1])),[2 1]);

save([DTfolderAll,'/gdat'],'gdat','-v7.3');

% %% Taking out noise at 100 & 200 
% [a,b] = butter( 4,[99/(srate/2), 101/(srate/2)],'stop');
% gdat = permute( filtfilt( a,b,permute( double( gdat ),[2 1])),[2 1]);
%  
% [a,b] = butter( 4,[199/(srate/2), 201/(srate/2)],'stop');
% gdat = permute( filtfilt( a,b,permute( double( gdat ),[2 1])),[2 1]);
% 
% save([DTfolderAll,'/gdat'],'gdat','-v7.3');

%% Convert to microvolts
if abs(gdat(1,1)) <0.001
    gdat = gdat*1e6;
end
save([DTfolderAll,'/gdat'],'gdat','-v7.3');

%% CAR.mat 
load([DTfolderAll,'/gdat']); 
Enum = 81;
%grouping = 16;
gdat_CAR = create_CAR(subj, block, Enum, bad_elecs, 'gdat', gdat);
% Note that elecs numbers will change again
save([DTfolderAll,'/gdat_CAR'],'gdat_CAR','-v7.3');
%% FFT plot for removing other noisy freqs
%sets up a structure called ecog
ecog.data = gdat_CAR;
ecog.srate = srate;

% here I included just the first channel which makes the estimation of the spectrum faster
ecog.data = ecog.data( 30,: );


% running the spectral analysis
tic;[S,F] = ecogSpectrum( ecog );toc;
x = squeeze( mean( S,3 ));

% plotting the results
plot( F,x )

%% Remove Bad Events

%define ictals 
ictals = [];


Events = markIctalEvents(ictals,Events,200,1000);
fprintf('Marked %i bad events...\n', sum([Events.badevent]));
fprintf('ReSaving Cleaned Events Structure...\n')
save([DTfolderAll,'/Events'],'Events');

%% Mark noisy elecs as zeros & save as gdat_CAR_good

gdat_CAR ([1:4 13:28 69],:)=0; %if we want to not look at some of the elecs
gdat_CAR_good=gdat_CAR;

save([DTfolderAll,'/gdat_CAR_good'],'gdat_CAR_good','-v7.3');

%% create new vectors of onsets for each condition without Bad events
onsetsViewR = [];
onsetsRestR = [];
onsetsGraspR = [];

for i = 1:size(Events,2)
    if 1 == strcmp(Events(i).cond,'viewR') && ~Events(i).badevent
        onsetsViewR = [onsetsViewR, Events(i).onset];
    elseif 1 == strcmp(Events(i).cond,'restR') && ~Events(i).badevent
        onsetsRestR = [onsetsRestR, Events(i).onset];
    elseif 1 == strcmp(Events(i).cond,'graspR') && ~Events(i).badevent
        onsetsGraspR = [onsetsGraspR, Events(i).onset];
    end;
end;
        
    
save([DTfolderAll,'/onsetsViewR'],'onsetsViewR');
save([DTfolderAll,'/onsetsRestR'],'onsetsRestR');
save([DTfolderAll,'/onsetsGraspR'],'onsetsGraspR');

%END OF PREPROCESSING BEFORE ANALYSIS

