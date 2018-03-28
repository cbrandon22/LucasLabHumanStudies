%Written by JFB
% Updates by AGR: 
% (1) automatically make good_leads
% (2) re-write jacksheet so it matches up with vox tool
% USER INPUTS
subj       = 'HUP128';
task       = 'corticalStim';
rawDirName = '';
badTags    = {'EKG','ekg','ECG','DC'}; % exclude from re-ref, Don't include num (e.g. use DC, not DC1)
%exclude these from jacksheet and re-ref
badLabels = {'Event','EKG1','EKG2','ekg1','ekg2','ECG1','ECG2','DC2','DC3','DC4','DC5',...
    'DC6','DC7','DC8','DC9','DC10','DC11','DC12','DC13','DC14','DC15',...
    'DC16','TRIG','OSAT','PR','Pleth','C8','C43','C44','C45','C46','C47','C48',...
    'C49','C50','C51','C52','C53','C54','C55','C56','C57','C58','C59',...
    'C60','C61','C62','C63','C64','C65','C66','C67','C68','C69','C70',...
    'C71','C72','C73','C74','C75','C76','C77','C78','C79','C80','C81',...
    'C82','C83','C84','C85','C86','C87','C88','C89','C90','C91','C92',...
    'C93','C94','C95','C96','C97','C98','C99','C100','C101','C102','C103',...
    'C104','C105','C106','C107','C108','C109','C110','C111','C112','C113',...
    'C114','C115','C116','C117','C118','C119','C120','C121','C122',...
    'C123','C124','C125','C126','C127','C128','C129','C130','C131','C132',...
    'C133','C134','C135','C136','C137','C138','C139','C140','C141','C142',...
    'C143','C144','C145','C146','C147','C148','C149','C150','C151','C152',...
    'C153','C154','C155','C156','C157','C158','C159','C160','C161','C162',...
    'C163','C164','C165','C166','C167','C168','C169','C170','C171','C172',...
    'C173','C174','C175','C176','C177','C178','C179','C180','C181','C182',...
    'C183','C184','C185','C186','C187','C188','C189','C190','C191','C192',...
    'C193','C194','C195','C196','C197','C198','C199','C200','C201','C202','C203',...
    'C204','C205','C206','C207','C208','C209','C210','C211','C212','C213',...
    'C214','C215','C216','C217','C218','C219','C220','C221','C222',...
    'C223','C224','C225','C226','C227','C228','C229','C230','C231','C232',...
    'C233','C234','C235','C236','C237','C238','C239','C240','C241','C242',...
    'C243','C244','C245','C246','C247','C248','C249','C250','C251','C252',...
    'C253','C254','C255','C256','C257','C258','C259','C260','C261','C262',...
    'C263','C264','C265','C266','C267','C268','C269','C270','C271','C272',...
    'C273','C274','C275','C276','C277','C278','C279','C280'}; 
% badLabels = {'PUSH','Event','PB','C43','C44','C45','C46','C47','C48',...
%     'C49','C50','C51','C52','C53','C54','C55','C56','C57','C58','C59',...
%     'C60','C61','C62','C63','C64','C65','C66','C67','C68','C69','C70',...
%     'C71','C72','C73','C74','C75','C76','C77','C78','C79','C80','C81',...
%     'C82','C83','C84','C85','C86','C87','C88','C89','C90','C91','C92',...
%     'C93','C94','C95','C96','C97','C98','C99','C100','C101','C102','C103'...
%     ,'C104','C105','C106','C107','C108','C109','C110','C111','C112','C113'...
%     ,'C114','C115','C116','C117','C118','C119','C120','C121','C122',...
%     'C123','C124','C125','C126','C127','C128','F7','T7','Tp11','P7','Fp1',...
%     'F3','C3','P3','O1','Fpz','Fz','Cz','Pz','Oz','Fp2','F4','C4','P4',...
%     'O2','F8','T8','Tpz','P8','TP11','TP12','DC1','DC2','DC3','DC4','DC5',...
%     'DC6','DC7','DC8','DC9','DC10','DC11','DC12','DC13','DC14','DC15',...
%     'DC16','TRIG','OSAT','PR','Pleth'}; 
%---------------------------------------
%---------------------------------------
%---------------------------------------
dirs = le_dirs(task);
dataDir  = dirs.data;
%dataDir = '/Volumes/Lucas_ECoG_Data/Consciousness/SleepData';
subjDir  = fullfile(dataDir,'eeg',subj);
%subjDir  = fullfile(dataDir,subj);
outDir   = fullfile(subjDir,'eeg.noreref');
rerefDir = fullfile(subjDir,'eeg.reref');
rawDir   = fullfile(subjDir,'raw',rawDirName);
docsDir  = fullfile(subjDir,'docs');
talDir   = fullfile(subjDir,'tal');
fprintf('\n')
eegFile  = getEDF_eeg_files(rawDir);
for k=1:length(eegFile)
  fprintf(' Splitting %d of %d files: \n\t%s',k,length(eegFile),eegFile{k})  
  edf_split(eegFile{k},subj,outDir,badLabels)
end
fprintf('\n\n')

% copy the jacksheet to the docs directory
jackMotherFile = fullfile(docsDir,'jacksheet.txt');
if ~exist(jackMotherFile,'file')
  fprintf('Making jacksheet in docs...')
  jackDaughter = fullfile(outDir,'jacksheet.txt');
  if ~exist(jackDaughter,'file');
    error(sprintf('\n\n\tNo jacksheet.txt in %s...not right\n\n',outDir))
  end
  system(['cp ' jackDaughter ' ' jackMotherFile]);
  pause(.1)
  fprintf('done\n')
else
  fprintf('Jacksheet exists in docs\n')
end

% make the electrode.m file
electrodesFile = fullfile(docsDir,'electrodes.m');
if ~exist(electrodesFile,'file')
  fprintf('Making electrodes.m in docs...')
  fid_jack = fopen(jackMotherFile,'r');
  X        = textscan(fid_jack,'%d%s');
  fclose(fid_jack);
  elecNum             = X{1};
  elecNam             = X{2};
  elecNam_stripped    = regexprep(elecNam,'\d+$','');
  un_elecNam_stripped = unique(elecNam_stripped);
  minMaxEls           = [];
  for k=1:length(un_elecNam_stripped)
    thisTag    = un_elecNam_stripped{k};
    indThisTag = strcmp(elecNam_stripped,thisTag);
    elsThisTag = elecNum(indThisTag);
    minMaxEls  = cat(1,minMaxEls,[min(elsThisTag) max(elsThisTag)]);
  end
  
  % now sort them
  [~,sortInd] = sort(minMaxEls(:,1));
  tagToWriteInElectrodesDotM = un_elecNam_stripped(sortInd);
  numToWriteInElectrodesDotM = minMaxEls(sortInd,:);
  
  % now write the electrodes.m file
  fid_elecs = fopen(electrodesFile,'w');
  %fprintf(fid_elecs,'%% electrodes.m file\n');
  %fprintf(fid_elecs,'%% Made within edf_split_wrapper.m for %s\n',subj);
  %fprintf(fid_elecs,'%%  ''''  using information from raw directory %s\n',rawDirName);
  %[~,user]=system('whoami');user=regexprep(user,'\n','');
  %fprintf(fid_elecs,'%%  ''''  by %s on %s\n',user,datestr(now));
  %fprintf(fid_elecs,'%%\n%%\n\n');
  fprintf(fid_elecs,'r=[\n');
  for k=1:size(numToWriteInElectrodesDotM,1)
    if all(cellfun(@isempty,regexp(tagToWriteInElectrodesDotM{k},badTags)))
    %if isempty(strfind(tagToWriteInElectrodesDotM{k},'EKG'));
      fprintf('  including %s in reference\n',tagToWriteInElectrodesDotM{k})
      fprintf(fid_elecs,'%3.0d , %3.0d;  %%%s  (OPTIONAL: fill in full name here)\n',...
	      numToWriteInElectrodesDotM(k,:),tagToWriteInElectrodesDotM{k});
    else
      fprintf('  **NOT** including %s in reference\n',tagToWriteInElectrodesDotM{k})
    end
  end
  fprintf(fid_elecs,'];\n\n');
  fclose(fid_elecs);  
  pause(.1)
  fprintf('done\n')

% no need to do this because edf_split has been updated - agr 6.29.2014
%   %Here, we will re-write the jacksheet from electrodes.m in the docs directory so that 
%   %it has electrode names that match up with vox tool
%   cd(docsDir)
%   !mv jacksheet.txt jacksheet_orig.txt
%   electrodes2jacksheet(subj)

 else
   fprintf('electrodes.m exists in docs\n')
 end


% make the leads.txt file
leadsFile = fullfile(talDir,'leads.txt');
if ~exist(leadsFile,'file')
  fprintf('Making leads.txt in docs...')
  fid_jack  = fopen(fullfile(docsDir,'jacksheet.txt'),'r');
  fid_leads = fopen(leadsFile,'w');
  X         = textscan(fid_jack,'%d%s');
  fclose(fid_jack);
  elecNum   = X{1};
  for e=1:length(elecNum)
    fprintf(fid_leads,'%d\n',elecNum(e));
  end  
  fclose(fid_leads);
  fprintf('done.\n')
  
  %write out a good_leads (can be edited later)
  cd(talDir)
  !cp leads.txt good_leads.txt

else
  fprintf('leads.txt exists in tal\n')
end

% reref it
fprintf('\nreref me: \n')
for k=1:length(eegFile)  
  thisFileExt = edf_split(eegFile{k},subj,outDir);
  tagNameFile = fullfile(docsDir,'electrodes.m');
  run(tagNameFile)
  reref({thisFileExt},r,rerefDir,talDir);
end

% mk the bipolar thing
%mkGridStructForBipolarFromJacksheet(subj,badTags)
