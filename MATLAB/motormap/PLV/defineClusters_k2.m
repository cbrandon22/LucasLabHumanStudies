% load k2_clusInfo
dirs = le_dirs('motormap');
clusInfoDir = fullfile(dirs.scratch,'PHASE/k2');
cd(clusInfoDir);
load('k2_clusInfo.mat');

clear analyze
useBipolar=0; %0 if using monopolar
allClus={'clus1' 'clus2' 'clus0' 'clusH'};
for analyze=1:5
    switch analyze
        case 1
            subj='HUP111';
            currVec=strcmp(subj,subjVec);
            currClusInfo=clus_info(currVec);
            count1=0;
            count2=0;
            count0=0;
            countH=0;
            for i=1:length(currClusInfo)
                if currClusInfo(i).cluster==1
                    disp('clus1')
                    count1=count1+1;
                    clus1chan(count1)=currClusInfo(i).eLbl1;
                    clus1chan(count1+1)=currClusInfo(i).eLbl2;
                    count1=count1+1;
                elseif currClusInfo(i).cluster==2
                    disp('clus2')
                    count2=count2+1;
                    clus2chan(count2)=currClusInfo(i).eLbl1;
                    clus2chan(count2+1)=currClusInfo(i).eLbl2;
                    count2=count2+1;
                elseif currClusInfo(i).cluster==0
                    disp('clus0')
                    count0=count0+1;
                    clus0chan(count0)=currClusInfo(i).eLbl1;
                    clus0chan(count0+1)=currClusInfo(i).eLbl2;
                    count0=count0+1;
                    if sum(strcmp(currClusInfo(i).anatAbbr,{'L-HIP','R-HIP','L-PHIP','R-PHIP'}))>0
                        disp('hippo')
                        countH=countH+1;
                        clusHchan(countH)=currClusInfo(i).eLbl1;
                        clusHchan(countH+1)=currClusInfo(i).eLbl2;
                        countH=countH+1;
                    end
                end
            end
            
        case 2
            subj='HUP084';
            currVec=strcmp(subj,subjVec);
            currClusInfo=clus_info(currVec);
            count1=0;
            count2=0;
            count0=0;
            countH=0;
            for i=1:length(currClusInfo)
                if currClusInfo(i).cluster==1
                    disp('clus1')
                    count1=count1+1;
                    clus1chan(count1)=currClusInfo(i).eLbl1;
                    clus1chan(count1+1)=currClusInfo(i).eLbl2;
                    count1=count1+1;
                elseif currClusInfo(i).cluster==2
                    disp('clus2')
                    count2=count2+1;
                    clus2chan(count2)=currClusInfo(i).eLbl1;
                    clus2chan(count2+1)=currClusInfo(i).eLbl2;
                    count2=count2+1;
                elseif currClusInfo(i).cluster==0
                    disp('clus0')
                    count0=count0+1;
                    clus0chan(count0)=currClusInfo(i).eLbl1;
                    clus0chan(count0+1)=currClusInfo(i).eLbl2;
                    count0=count0+1;
                    if sum(strcmp(currClusInfo(i).anatAbbr,{'L-HIP','R-HIP','L-PHIP','R-PHIP'}))>0
                        disp('hippo')
                        countH=countH+1;
                        clusHchan(countH)=currClusInfo(i).eLbl1;
                        clusHchan(countH+1)=currClusInfo(i).eLbl2;
                        countH=countH+1;
                    end
                end
            end
            
        case 3
            subj='HUP087';
            currVec=strcmp(subj,subjVec);
            currClusInfo=clus_info(currVec);
            count1=0;
            count2=0;
            count0=0;
            countH=0;
            for i=1:length(currClusInfo)
                if currClusInfo(i).cluster==1
                    disp('clus1')
                    count1=count1+1;
                    clus1chan(count1)=currClusInfo(i).eLbl1;
                    clus1chan(count1+1)=currClusInfo(i).eLbl2;
                    count1=count1+1;
                elseif currClusInfo(i).cluster==2
                    disp('clus2')
                    count2=count2+1;
                    clus2chan(count2)=currClusInfo(i).eLbl1;
                    clus2chan(count2+1)=currClusInfo(i).eLbl2;
                    count2=count2+1;
                elseif currClusInfo(i).cluster==0
                    disp('clus0')
                    count0=count0+1;
                    clus0chan(count0)=currClusInfo(i).eLbl1;
                    clus0chan(count0+1)=currClusInfo(i).eLbl2;
                    count0=count0+1;
                    if sum(strcmp(currClusInfo(i).anatAbbr,{'L-HIP','R-HIP','L-PHIP','R-PHIP'}))>0
                        disp('hippo')
                        countH=countH+1;
                        clusHchan(countH)=currClusInfo(i).eLbl1;
                        clusHchan(countH+1)=currClusInfo(i).eLbl2;
                        countH=countH+1;
                    end
                end
            end
            
        case 4
            subj='HUP091';
            currVec=strcmp(subj,subjVec);
            currClusInfo=clus_info(currVec);
            count1=0;
            count2=0;
            count0=0;
            countH=0;
            for i=1:length(currClusInfo)
                if currClusInfo(i).cluster==1
                    disp('clus1')
                    count1=count1+1;
                    clus1chan(count1)=currClusInfo(i).eLbl1;
                    clus1chan(count1+1)=currClusInfo(i).eLbl2;
                    count1=count1+1;
                elseif currClusInfo(i).cluster==2
                    disp('clus2')
                    count2=count2+1;
                    clus2chan(count2)=currClusInfo(i).eLbl1;
                    clus2chan(count2+1)=currClusInfo(i).eLbl2;
                    count2=count2+1;
                elseif currClusInfo(i).cluster==0
                    disp('clus0')
                    count0=count0+1;
                    clus0chan(count0)=currClusInfo(i).eLbl1;
                    clus0chan(count0+1)=currClusInfo(i).eLbl2;
                    count0=count0+1;
                    if sum(strcmp(currClusInfo(i).anatAbbr,{'L-HIP','R-HIP','L-PHIP','R-PHIP'}))>0
                        disp('hippo')
                        countH=countH+1;
                        clusHchan(countH)=currClusInfo(i).eLbl1;
                        clusHchan(countH+1)=currClusInfo(i).eLbl2;
                        countH=countH+1;
                    end
                end
            end
            
        case 5
            subj='HUP092';
            currVec=strcmp(subj,subjVec);
            currClusInfo=clus_info(currVec);
            count1=0;
            count2=0;
            count0=0;
            countH=0;
            for i=1:length(currClusInfo)
                if currClusInfo(i).cluster==1
                    disp('clus1')
                    count1=count1+1;
                    clus1chan(count1)=currClusInfo(i).eLbl1;
                    clus1chan(count1+1)=currClusInfo(i).eLbl2;
                    count1=count1+1;
                elseif currClusInfo(i).cluster==2
                    disp('clus2')
                    count2=count2+1;
                    clus2chan(count2)=currClusInfo(i).eLbl1;
                    clus2chan(count2+1)=currClusInfo(i).eLbl2;
                    count2=count2+1;
                elseif currClusInfo(i).cluster==0
                    disp('clus0')
                    count0=count0+1;
                    clus0chan(count0)=currClusInfo(i).eLbl1;
                    clus0chan(count0+1)=currClusInfo(i).eLbl2;
                    count0=count0+1;
                    if sum(strcmp(currClusInfo(i).anatAbbr,{'L-HIP','R-HIP','L-PHIP','R-PHIP','L-AMYG', 'R-AMYG'}))>0
                        disp('hippo')
                        countH=countH+1;
                        clusHchan(countH)=currClusInfo(i).eLbl1;
                        clusHchan(countH+1)=currClusInfo(i).eLbl2;
                        countH=countH+1;
                    end
                end
            end
            
    end
    run vb_clusEEG.m
end








