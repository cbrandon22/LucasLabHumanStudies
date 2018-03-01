function gdat_CAR = create_CAR(subj,Enum,bad_elecs, varargin)
%   create_CAR(subj,block,Enum,grouping,'gropuLUT', groupLUT, 'gdat',gdat, 'exclude',exclude)
%     The function will create a CAR data matrix and save it.
%       subj    - subj code
%       block   - block name
%       Enum    - Number of electrodes in the matrix
%       grouping- Size of groups to average
%       exclude - More electrodes to exclude (adds them to bad_elecs group)
%       groupLUT- Manually enter the electrode grouping.
%       gdat - provide gdat instead of loading it in from file
%   The function will first group create a CAR in groups of grouping and
%   then a common CAR for all valid electrodes (in elecs). The final CAR
%

   % get_subj_globals(subj,block);
   
    % Set defaults
    for n=1:2:length(varargin)-1
        switch lower(varargin{n})
            case 'exclude'
                exclude = varargin{n+1};
                bad_elecs = union(bad_elecs,exclude);
                fprintf('Excluding %s\n',num2str(bad_elecs));
            case 'gdat'
                gdat = varargin{n+1};
            case 'grouplut'
                groupLUT = varargin{n+1};
            case 'grouping'
                grouping = varargin{n+1};
        end
    end

    if ~(exist('gdat','var'))
        load 'gdat';
    end

    if ~(exist('grouping','var'))
        grouping = Enum;
    end

    chRank = ones(1,Enum);
    chRank(bad_elecs) = 0;

    Ngroups = ceil(Enum/grouping)
    % Ngroups = ceil(elecs(end)/grouping);
    %Ngroups = ceil(length(chRank)/grouping); %edit mh

    disp('CAR started...');
    if (~exist('groupLUT'))
        groups=[];
        for cnt = 1:Ngroups        % create grouping lookup table
            groups = [groups cnt*ones(1,grouping)];
        end
        %Ngroups = ceil(elecs(end)/grouping);
        %Ngroups = ceil(length(chRank)/grouping); %edit mh
    else
        groups = groupLUT;
        Ngroups = max(groupLUT);
    end
    CAR = zeros(Ngroups,length(gdat));
    CAR_all = zeros(1,length(gdat));
    disp('Creating group CAR');

    for e = 1:numel(groups) %changed from elecs
        gdat(e,:) = gdat(e,:) - mean(gdat(e,:));    % remove mean before CAR (only for good elecs)

        if chRank(e)
            CAR(groups(e),:) = CAR(groups(e),:) + gdat(e,:);        % sum CAR in groups (banks)
        end
    end

    %for cnt = 1:Ngroups
    %        if grouping*cnt>length(chRank)
    %            CAR(cnt,:) = CAR(cnt,:)./sum(chRank(grouping*(cnt-1)+1:end));  %divide by number of valid channels in the group
            %else
            %    CAR(cnt,:) = CAR(cnt,:)./sum(chRank(grouping*(cnt-1)+1:grouping*cnt));  %divide by number of valid channels in the group
                %easier to do find(groups == cnt) %mh
            %end
    %end

    %edit MH to allow for groups/banks of different sizes - based on lookup table. makes grouping argument irrelevant.
    for cnt = 1:Ngroups
        CAR(cnt,:) = CAR(cnt,:)./sum(chRank(groups == cnt));
    end

    disp('Removing group CAR, elecs');
    for e = 1:Enum
        gdat(e,:) = gdat(e,:) - CAR(groups(e),:);             % remove group CAR groups(e)?
        CAR_all = CAR_all + gdat(e,:);                        % sum CAR for all channels
    end
    disp('Removing group CAR, bad_elecs');
    for e = bad_elecs
        gdat(e,:) = gdat(e,:) - CAR(groups(e),:);             % remove group CAR for bad_elecs
    end
    CAR_all = CAR_all./Enum;

    %disp('saving CAR, CAR_group');
    %save([pwd,'\SFC1\CAR'],'CAR*','chRank');
    %save([pwd,'\SFC1\gdat_CAR_group'],'gdat','-v7.3');
    %disp('Removing total CAR');
   % for e = [elecs bad_elecs]
    %    gdat(e,:) = gdat(e,:) - CAR_all;                      % remove CAR for all elecs
    %end
    %disp('saving gdat_CAR');
    gdat_CAR = gdat;
    %save([pwd,'\SFC1\gdat_CAR'],'gdat','-v7.3');
