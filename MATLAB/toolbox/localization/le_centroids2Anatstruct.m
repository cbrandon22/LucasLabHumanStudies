function[anatStruct_bipol,anatStruct_mono] = le_centroids2Anatstruct(subj)
% This function loads data from the SUBJ_labels.csv file and generates an
% anatStruct that can be used for analysis
% Inputs
% subj = 'HUP001';
% bipolarFlag = 1;
task = 'motormap';


% get dirs
d = le_dirs(task);
count = 0;

% open jacksheet
%jacFile = fullfile(d.data,'eeg',task,subj,'docs','jacksheet.txt');
jacFile = fullfile(d.data,'eeg',subj,'docs','jacksheet.txt');
fid = fopen(jacFile,'r');
JAC=textscan(fid,'%d%s','delimiter','\t');
JAC{:,2} = strtrim(JAC{:,2});
fclose(fid);

% load labels file
%fname = fullfile(d.data,'eeg',task,subj,'tal',[subj '_labels.csv']);
fname = fullfile(d.data,'eeg',subj,'tal',[subj '_labels.csv']);
[fid] = fopen(fname,'r');

% load anatomical number label pairings
anatFile = fullfile(d.reportDir,'anatLbls.csv');
[anatid] = fopen(anatFile,'r');

% get anatLine
anatLine = fgetl(anatid);
lineNum = 1;

% build anatCell
anatCell = cell(120,3);

while ischar(anatLine)
    % extract line into cell array
    A = textscan(anatLine,'%f%s%s','delimiter',',','EmptyValue',nan);
    
    % add line to anatCell
    anatCell(lineNum,:) = {A{1} A{2} A{3}};
    
    % next line
    anatLine = fgetl(anatid);
    lineNum = lineNum+1;
end
fclose(anatid);

%read data file
count = 0;anatStruct_bipol = struct;anatStruct_mono = struct;

while true
    % get line, break if empty(end of file)
    tline = fgetl(fid);
    if ~ischar(tline)
        break
    end
    
    % extract line into cell array
    C = textscan(tline,...
    '%f%f%f%s%f%f%f','delimiter',',','EmptyValue',nan);

    % get lbl
    elecLbl = C{4};
    anatNum = C{7};
    anatLbl = anatCell{anatNum,2};
    anatAbbr = anatCell{anatNum,3};
    
    % get elecnum using jacksheet
    idx = strcmp(JAC{:,2},elecLbl);
    
    if sum(idx)>0
    % update counter
        count = count + 1; 
        
        %update struct
        anatStruct_mono(count).elecNum = JAC{1}(idx);
        anatStruct_mono(count).elecLbl = elecLbl{1};
        anatStruct_mono(count).eLbl = num2str(anatStruct_mono(count).elecNum);
        anatStruct_mono(count).X = C{1};
        anatStruct_mono(count).Y = C{2};
        anatStruct_mono(count).Z = C{3};
        anatStruct_mono(count).XYZ = [C{1} C{2} C{3}];
        anatStruct_mono(count).dim = [C{5} C{6}];
        anatStruct_mono(count).anatNum = anatNum;
        anatStruct_mono(count).anatLbl = anatLbl{1};
        anatStruct_mono(count).anatAbbr = anatAbbr{1};
    end
end

%create bipolar structure by finding all possible pairs of monopolar structure
count = 0;
for x = 1:length(anatStruct_mono)
    for y = 1:length(anatStruct_mono)
        
        
        if x >= y
            continue
        end
        
        count = count + 1;
        anatStruct_bipol(count).elecNum = [anatStruct_mono(x).elecNum anatStruct_mono(y).elecNum];
        anatStruct_bipol(count).elecLbl = [anatStruct_mono(x).elecLbl '-' anatStruct_mono(y).elecLbl];
        anatStruct_bipol(count).eLbl = [anatStruct_mono(x).eLbl '-' anatStruct_mono(y).eLbl];
        anatStruct_bipol(count).e1_idx = x;
        anatStruct_bipol(count).e2_idx = y;
        anatStruct_bipol(count).XYZ = bipolarHalfDistance_local(anatStruct_mono(x).XYZ,...
            anatStruct_mono(y).XYZ);
        anatStruct_bipol(count).X = anatStruct_bipol(count).XYZ(1);
        anatStruct_bipol(count).Y = anatStruct_bipol(count).XYZ(2);
        anatStruct_bipol(count).Z = anatStruct_bipol(count).XYZ(3);
        anatStruct_bipol(count).distance = sqrt(sum((anatStruct_mono(x).XYZ-anatStruct_mono(y).XYZ).^2));
        
    end
end

% filter bipolar pairs that have a distance between 8 and 12 mm
idx_to_keep = [anatStruct_bipol.distance]<=10.5 & [anatStruct_bipol.distance]>=9.5;
anatStruct_bipol = anatStruct_bipol(idx_to_keep);

% if bpCentroid_labels have been generated, add anatomical info
bpAnatFile = fullfile(d.data,'eeg',subj,'tal',[subj '_bpCentroid_labels.csv']);
lineNum = 1;
if exist(bpAnatFile)
    % load bipolar labels
    [bpanat] = fopen(bpAnatFile,'r');
    bpAnatLine = fgetl(bpanat);
    % add anatomical info to anatStruct_bipol
    while ischar(bpAnatLine)
    % extract line into cell array
    C = textscan(bpAnatLine,'%f%f%f%f','delimiter',',','EmptyValue',nan);
    
    % get anat info
    anatNum = C{4};
    anatLbl = anatCell{anatNum,2};
    anatAbbr = anatCell{anatNum,3};
    anatStruct_bipol(lineNum).anatNum = anatNum;
    anatStruct_bipol(lineNum).anatLbl = anatLbl{1};
    anatStruct_bipol(lineNum).anatAbbr = anatAbbr{1};
    % add anat info to bipol struct
    

    % next line
    bpAnatLine = fgetl(bpanat);
    lineNum = lineNum+1;
    end
else
    warning('No bipolar anatomical lables')
end

%keyboard

function xyz_out = bipolarHalfDistance_local(XYZ1,XYZ2)
tmp = -XYZ1+XYZ2;
halftmp = tmp./2;
xyz_out = XYZ1 + halftmp;
