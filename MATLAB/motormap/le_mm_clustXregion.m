function le_mm_clustXregion
% function clustXregion
%   Hypothesis test of motor task data.
%   Ho: ECoG power profile clusters are uniformally distributed across
%   brain regions.
%
%   DR 07/2016

% parameters
dirs = le_dirs('motormap');
ddir = dirs.reportDir;
% regL = {'L-PRE','R-PRE','L-POST','R-POST',[],[]; % region 1 label(s)
%         'L-T1','R-T1','L-T2','R-T2','L-T3','R-T3'}; % region 2 label(s)
%                                                     % region 3 label(s) ...
% regL = {'L-PRE','R-PRE',[],[],[],[],[],[],[],[];
%         'L-F1','R-F1',[],[],[],[],[],[],[],[];
%         'L-F2','R-F2',[],[],[],[],[],[],[],[];
%         'L-F3OP','R-F3OP','L-F3T','R-F3T',[],[],[],[],[],[];
%         'L-RO','R-RO',[],[],[],[],[],[],[],[];
%         'L-GR','R-GR',[],[],[],[],[],[],[],[];
%         'L-HIPP','R-HIPP',[],[],[],[],[],[],[],[];
%         'L-O2','R-O2','L-O3','R-O3',[],[],[],[],[],[];
%         'L-FUSI','R-FUSI',[],[],[],[],[],[],[],[];
%         'L-POST','R-POST',[],[],[],[],[],[],[],[];
%         'L-P1','R-P1','L-P2','R-P2',[],[],[],[],[],[];
%         'L-SMG','R-SMG',[],[],[],[],[],[],[],[];
%         'L-AG','R-AG',[],[],[],[],[],[],[],[];
%         'L-T1','R-T1','L-T1P','R-T1P','L-T2','R-T2','L-T2P','R-T2P','L-T3','R-T3'};
% regL = {'L-PRE','R-PRE',[],[],[],[],[],[],[],[];
%         'L-POST','R-POST',[],[],[],[],[],[],[],[];
%         'L-F2','R-F2',[],[],[],[],[],[],[],[];
%         'L-T1','R-T1','L-T1P','R-T1P','L-T2','R-T2','L-T2P','R-T2P','L-T3','R-T3'};
%         
% regL = {'L-PRE','R-PRE','L-F1','R-F1','L-F2','R-F2','L-F3OP','R-F3OP','L-F3T','R-F3T','L-RO','R-RO','L-GR','R-GR';
%         'L-HIPP','R-HIPP','L-FUSI','R-FUSI','L-T1','R-T1','L-T1P','R-T1P','L-T2','R-T2','L-T2P','R-T2P','L-T3','R-T3';
%         'L-O2','R-O2','L-O3','R-O3',[],[],[],[],[],[],[],[],[],[];
%         'L-POST','R-POST','L-P1','R-P1','L-P2','R-P2','L-SMG','R-SMG','L-AG','R-AG',[],[],[],[]};
% regL = {'L-T1','L-T1P','L-T2','L-T2P','L-T3';
%         'R-T1','R-T1P','R-T2','R-T2P','R-T3'};

% load spreadsheet data
cd(ddir);
%[~,labl] = xlsread('anatClusters.xls','A1:A1000');
%clust = xlsread('anatClusters.xls','B1:B1000');
%%clust = xlsread('anatSubclustersClus3.xls','B1:B1000');
%[~,labl] = xlsread('anatSubclustersClus1.xls','A1:A1000');
%clust = xlsread('anatSubclustersClus1.xls','B1:B1000');
[~,labl] = xlsread('anat_clusterk22.xlsx','A1:A1000');%194 is end of motor bps
clust = xlsread('anat_clusterk22.xlsx','B1:B1000');%194 is end of motor bps

regL = unique(labl)';
motorReg = {'L-PRE','R-PRE','L-POST','R-POST'};
Lia = ismember(regL,motorReg);
regL = regL(~Lia);
regL(2,1:length(motorReg)) = motorReg;

% region categories
N = size(regL);
region = zeros(length(clust),N(1));
for ii = 1:N(1)
    for jj = 1:N(2)
        region(:,ii) = region(:,ii) | cellfun(@(x) ~isempty(x),strfind(labl,regL{ii,jj}));
    end
    region(:,ii) = region(:,ii)*ii; %comment out to remove regions with few bps
    % Remove regions with fewer than 15 significant electrodes
%     if sum(region(:,ii)) >= 15
%         region(:,ii) = region(:,ii)*ii;
%     else
%         region(:,ii) = 0;
%     end
end

region = sum(region,2);
clust(region==0) = []; % remove data not in any of the specified regions
region(region==0) = [];

% contingency table and chi-squared test
[table,chi2,p] = crosstab(clust,region)
