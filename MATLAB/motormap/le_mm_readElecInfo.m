function [eInfoStruct] = le_mm_readElecInfo
% Read ElectrodeInfo.xls and bring in variables
dirs = le_dirs;
cd(dirs.scratch);
[~,subjects] = xlsread('ElectrodeInfo.xlsx','A1:A1000');
[~,elecLbl1] = xlsread('ElectrodeInfo.xlsx','B1:B1000');
[~,elecLbl2] = xlsread('ElectrodeInfo.xlsx','C1:C1000');
eLbl1 = xlsread('ElectrodeInfo.xlsx','D1:D1000');
eLbl2 = xlsread('ElectrodeInfo.xlsx','E1:E1000');
X = xlsread('ElectrodeInfo.xlsx','F1:F1000');
Y = xlsread('ElectrodeInfo.xlsx','G1:G1000');
Z = xlsread('ElectrodeInfo.xlsx','H1:H1000');
[~,anatAbbr] = xlsread('ElectrodeInfo.xlsx','I1:I1000');
cluster = xlsread('ElectrodeInfo.xlsx','J1:J1000');
for i = length(subjects):-1:1
    eInfoStruct(i).subj = subjects{i};
    eInfoStruct(i).elecLbl1 = elecLbl1{i};
    eInfoStruct(i).elecLbl2 = elecLbl2{i};
    eInfoStruct(i).eLbl1 = eLbl1(i);
    eInfoStruct(i).eLbl2 = eLbl2(i);
    eInfoStruct(i).X = X(i);
    eInfoStruct(i).Y = Y(i);
    eInfoStruct(i).Z = Z(i);
    eInfoStruct(i).anatAbbr = anatAbbr{i};
    eInfoStruct(i).cluster = cluster(i);
end
end
