hold on
xlim([-2 6.25])
ylim([-1.5 1.5])
plot([0 0], [-7 7], 'b-', 'linewidth', 2)
plot([-3 7], [0 0], 'k--', 'linewidth', 2)
plot([0 5], [-1 -1], 'k-', 'linewidth' ,2)
xlabel('Time (s)')
ylabel('Change in PLV')
legend Move Wait



figure
pcolor(squeeze(mean(HUP087_clus1_plvGammaWait,1)))
shading flat; colorbar
caxis ([0 1])

mean2(HUP087_clus2_plvThetaMove(1025:end,:,:)) %average the PLV from after go-cue across all regions
mean2(HUP087_clus2_plvThetaWait(1025:end,:,:)) 
mean2(HUP084_clus1_plvGammaMove(1025:end,:,:)) 
mean2(HUP084_clus1_plvGammaWait(1025:end,:,:)) 

%all clus 1 ROIs
allSubj_clus1_normGammaMove(1,:)=HUP087_clus1_m2_avg_normGammaMove_LF2_LT3; %9
allSubj_clus1_normGammaMove(2,:)=HUP087_clus1_m2_avg_normGammaMove_LPOST_LF2; %9
allSubj_clus1_normGammaMove(3,:)=HUP087_clus1_m2_avg_normGammaMove_LPOST_LT3; %9
allSubj_clus1_normGammaMove(4,:)=HUP087_clus1_m2_avg_normGammaMove_LPRE_LF2; %9
allSubj_clus1_normGammaMove(5,:)=HUP087_clus1_m2_avg_normGammaMove_LPRE_LPOST; %9
allSubj_clus1_normGammaMove(6,:)=HUP087_clus1_m2_avg_normGammaMove_LPRE_LT3; %9
allSubj_clus1_normGammaMove(7,:)=HUP090_clus1_m2_avg_normGammaMove_LPRE_LF2; %1
allSubj_clus1_normGammaMove(8,:)=HUP091_clus1_m2_avg_normGammaMove_LF2_LT3; %3
allSubj_clus1_normGammaMove(9,:)=HUP091_clus1_m2_avg_normGammaMove_LPOST_LF2; %3
allSubj_clus1_normGammaMove(10,:)=HUP091_clus1_m2_avg_normGammaMove_LPOST_LT3; %3
allSubj_clus1_normGammaMove(11,:)=HUP092_clus1_m2_avg_normGammaMove_LF2_LT3; %6
allSubj_clus1_normGammaMove(12,:)=HUP092_clus1_m2_avg_normGammaMove_LPOST_LF2; %18
allSubj_clus1_normGammaMove(13,:)=HUP092_clus1_m2_avg_normGammaMove_LPOST_LT3; %3
allSubj_clus1_normGammaMove(14,:)=HUP092_clus1_m2_avg_normGammaMove_LPRE_LF2; %9
allSubj_clus1_normGammaMove(15,:)=HUP092_clus1_m2_avg_normGammaMove_LPRE_LPOST; %3
allSubj_clus1_normGammaMove(16,:)=HUP092_clus1_m2_avg_normGammaMove_LPRE_LT3; %3  TOTAL pairwise n=106, across 4 subjects



allSubj_clus1_normThetaMove_motor(1,:)=HUP087_clus1_m2_avg_normThetaMove_LPOST_LF2; %9
allSubj_clus1_normThetaMove_motor(2,:)=HUP087_clus1_m2_avg_normThetaMove_LPRE_LF2; %9
allSubj_clus1_normThetaMove_motor(3,:)=HUP087_clus1_m2_avg_normThetaMove_LPRE_LPOST; %9
allSubj_clus1_normThetaMove_motor(4,:)=HUP090_clus1_m2_avg_normThetaMove_LPRE_LF2; %1
allSubj_clus1_normThetaMove_motor(5,:)=HUP091_clus1_m2_avg_normThetaMove_LPOST_LF2; %3
allSubj_clus1_normThetaMove_motor(6,:)=HUP092_clus1_m2_avg_normThetaMove_LPOST_LF2; %18
allSubj_clus1_normThetaMove_motor(7,:)=HUP092_clus1_m2_avg_normThetaMove_LPRE_LF2; %18
allSubj_clus1_normThetaMove_motor(8,:)=HUP092_clus1_m2_avg_normThetaMove_LPRE_LPOST; %9 TOTAL pairwise n=76, across 4 subjects






%all clus 2 ROIs, happen to be only in motor network
allSubj_clus2_normGammaWait(1,:)=HUP087_clus2_m2_avg_normGammaWait_LPRE_LPOST; %7
allSubj_clus2_normGammaWait(2,:)=HUP091_clus2_m2_avg_normGammaWait_LPRE_LPOST; %15   TOTAL  pairwise n=22, across 2 subjects




%all clus3
allSubj_clus3_normThetaWait(1,:)=HUP084_clus3_m2_avg_normThetaWait_LPOST_LT3; %pairwise n=1
allSubj_clus3_normThetaWait(2,:)=HUP087_clus3_m2_avg_normThetaWait_LPRE_LPOST; %45
allSubj_clus3_normThetaWait(3,:)=HUP089_clus3_m2_avg_normThetaWait_RPOST_RF2; %1
allSubj_clus3_normThetaWait(4,:)=HUP089_clus3_m2_avg_normThetaWait_RPRE_RF2; %4
allSubj_clus3_normThetaWait(5,:)=HUP089_clus3_m2_avg_normThetaWait_RPRE_RPOST; %4
allSubj_clus3_normThetaWait(6,:)=HUP091_clus3_m2_avg_normThetaWait_LF2_LT3; %8
allSubj_clus3_normThetaWait(7,:)=HUP091_clus3_m2_avg_normThetaWait_LPOST_LF2; %16
allSubj_clus3_normThetaWait(8,:)=HUP091_clus3_m2_avg_normThetaWait_LPOST_LT3; %8
allSubj_clus3_normThetaWait(9,:)=HUP091_clus3_m2_avg_normThetaWait_LPRE_LF2; %12
allSubj_clus3_normThetaWait(10,:)=HUP091_clus3_m2_avg_normThetaWait_LPRE_LPOST; %12
allSubj_clus3_normThetaWait(11,:)=HUP091_clus3_m2_avg_normThetaWait_LPRE_LT3; %6
allSubj_clus3_normThetaWait(12,:)=HUP092_clus3_m2_avg_normThetaWait_LPOST_LF2; %6
allSubj_clus3_normThetaWait(13,:)=HUP092_clus3_m2_avg_normThetaWait_LPRE_LF2; %5
allSubj_clus3_normThetaWait(14,:)=HUP092_clus3_m2_avg_normThetaWait_LPRE_LPOST; %30
allSubj_clus3_normThetaWait(15,:)=HUP111_clus3_m2_avg_normThetaWait_RPRE_RPOST; %40   TOTAL pairwise n=198, across 6 subj



allSubj_clus3_normGammaWait_motor(1,:)=HUP087_clus3_m2_avg_normGammaWait_LPRE_LPOST; %45
allSubj_clus3_normGammaWait_motor(2,:)=HUP089_clus3_m2_avg_normGammaWait_RPOST_RF2; %1
allSubj_clus3_normGammaWait_motor(3,:)=HUP089_clus3_m2_avg_normGammaWait_RPRE_RF2; %4
allSubj_clus3_normGammaWait_motor(4,:)=HUP089_clus3_m2_avg_normGammaWait_RPRE_RPOST; %4
allSubj_clus3_normGammaWait_motor(5,:)=HUP091_clus3_m2_avg_normGammaWait_LPOST_LF2; %16
allSubj_clus3_normGammaWait_motor(6,:)=HUP091_clus3_m2_avg_normGammaWait_LPRE_LF2; %12
allSubj_clus3_normGammaWait_motor(7,:)=HUP091_clus3_m2_avg_normGammaWait_LPRE_LPOST; %12
allSubj_clus3_normGammaWait_motor(8,:)=HUP092_clus3_m2_avg_normGammaWait_LPOST_LF2; %6
allSubj_clus3_normGammaWait_motor(9,:)=HUP092_clus3_m2_avg_normGammaWait_LPRE_LF2; %5
allSubj_clus3_normGammaWait_motor(10,:)=HUP092_clus3_m2_avg_normGammaWait_LPRE_LPOST; %30
allSubj_clus3_normGammaWait_motor(11,:)=HUP111_clus3_m2_avg_normGammaWait_RPRE_RPOST; %40    TOTAL pairwise n=175, across 5 subjects



%means

allSubj_mean_clus3_normThetaMove=mean(allSubj_clus3_normThetaMove);
allSubj_mean_clus3_normThetaWait=mean(allSubj_clus3_normThetaWait);
allSubj_mean_clus3_normGammaWait=mean(allSubj_clus3_normGammaWait);
allSubj_mean_clus3_normGammaMove=mean(allSubj_clus3_normGammaMove);


allSubj_std_clus3_normThetaMove=std(allSubj_clus3_normThetaMove);
allSubj_std_clus3_normThetaWait=std(allSubj_clus3_normThetaWait);
allSubj_std_clus3_normGammaWait=std(allSubj_clus3_normGammaWait);
allSubj_std_clus3_normGammaMove=std(allSubj_clus3_normGammaMove);



allSubj_mean_clus1_normThetaMove_motor=mean(allSubj_clus1_normThetaMove_motor);
allSubj_mean_clus1_normThetaWait_motor=mean(allSubj_clus1_normThetaWait_motor);
allSubj_mean_clus1_normGammaWait_motor=mean(allSubj_clus1_normGammaWait_motor);
allSubj_mean_clus1_normGammaMove_motor=mean(allSubj_clus1_normGammaMove_motor);


allSubj_std_clus1_normThetaMove_motor=std(allSubj_clus1_normThetaMove_motor);
allSubj_std_clus1_normThetaWait_motor=std(allSubj_clus1_normThetaWait_motor);
allSubj_std_clus1_normGammaWait_motor=std(allSubj_clus1_normGammaWait_motor);
allSubj_std_clus1_normGammaMove_motor=std(allSubj_clus1_normGammaMove_motor);




%% all pairwise

% count=0; %only for first subj, then comment out
   clear i j
    for i=1:size(HUP084_clus2_m2_normThetaMove,2)
        for j=1:size(HUP084_clus2_m2_normThetaMove,3)
            if i<j
                count=count+1;
                blahTM(count,:)=HUP084_clus2_m2_normThetaMove(:,i,j)';
                blahTW(count,:)=HUP084_clus2_m2_normThetaWait(:,i,j)';
                blahGM(count,:)=HUP084_clus2_m2_normGammaMove(:,i,j)';
                blahGW(count,:)=HUP084_clus2_m2_normGammaWait(:,i,j)';                
            else
                continue
            end
        end
    end
    clear i j
    

    
allSubj_PW_clus3_m2_normThetaMove=blahTM;
allSubj_PW_clus3_m2_normThetaWait=blahTW;
allSubj_PW_clus3_m2_normGammaMove=blahGM;
allSubj_PW_clus3_m2_normGammaWait=blahGW;

allSubj_PW_clus3_m2_avg_normThetaMove=mean(allSubj_PW_clus3_m2_normThetaMove);
allSubj_PW_clus3_m2_avg_normThetaWait=mean(allSubj_PW_clus3_m2_normThetaWait);
allSubj_PW_clus3_m2_avg_normGammaMove=mean(allSubj_PW_clus3_m2_normGammaMove);
allSubj_PW_clus3_m2_avg_normGammaWait=mean(allSubj_PW_clus3_m2_normGammaWait);

allSubj_PW_clus3_m2_std_normThetaMove=std(allSubj_PW_clus3_m2_normThetaMove);
allSubj_PW_clus3_m2_std_normThetaWait=std(allSubj_PW_clus3_m2_normThetaWait);
allSubj_PW_clus3_m2_std_normGammaMove=std(allSubj_PW_clus3_m2_normGammaMove);
allSubj_PW_clus3_m2_std_normGammaWait=std(allSubj_PW_clus3_m2_normGammaWait);

clear blah* count



%% significance testing 
%need vector of p-values created by movingwindow 190 samples, with 50
%sample overlap (approx .4 seconds with .1 second overlap)

%vs. baseline first to see if signal is modulated by stimulus
clear time overlap h1 h2 p1 p2 hob1 pob1 hob2 pob2
overlap=0;
bin=25; %samples
basetime=floor(475:bin:1525);
exptime=floor(2000:bin:3050);
sig1=mean(sigTM);
sig2=mean(sigTW);
alpha=0.05;


[hob1,pob1]=ttest(sig1(basetime(1):bin:basetime(end)),sig1(exptime(1):bin:exptime(end)), 'alpha', alpha);
[hob2,pob2]=ttest(sig2(basetime(1):bin:basetime(end)),sig2(exptime(1):bin:exptime(end)), 'alpha', alpha);

for i=1:length(basetime)-1
    if i==1;
        [h1(i),p1(i)]=ttest(sig1(basetime(i):basetime(i+1)),sig1(exptime(i):exptime(i+1)), 'alpha', alpha);
        [h2(i),p1(i)]=ttest(sig2(basetime(i):basetime(i+1)),sig2(exptime(i):exptime(i+1)), 'alpha', alpha);
        
    else
        [h1(i),p1(i)]=ttest(sig1(basetime(i)-overlap:basetime(i+1)-overlap),sig1(exptime(i)-overlap:exptime(i+1)-overlap), 'alpha', alpha);
        [h2(i),p2(i)]=ttest(sig2(basetime(i)-overlap:basetime(i+1)-overlap),sig2(exptime(i)-overlap:exptime(i+1)-overlap), 'alpha', alpha);
        
    end
end    

clear time overlap h p ho po
if hob1 && hob2
    % if passes baseline test, move on to this
    overlap=50;
    time=floor(1550:150:4620);
    alpha=0.01;
    for i=1:length(time)-1
        if i==1;
            [h(i),p(i)]=ttest(sig1(time(i):time(i+1)),sig2(time(i):time(i+1)), 'alpha', alpha);
        else
            [h(i),p(i)]=ttest(sig1(time(i)-overlap:time(i+1)-overlap),sig2(time(i)-overlap:time(i+1)-overlap), 'alpha', alpha);
        end
        
    end
    [ho,po]=ttest(sig1(time(1):time(end)),sig2(time(1):time(end)), 'alpha', alpha);
    if ho
        disp('*** Significant stimulus-induced connectivity modulation between both signals')
    end
else
    disp('No significant stimulus-induced modulation of one or both signals')
    disp(pob1)
    disp(pob2)
end


%% extract batchLoad
xxStruct=allSubj_k2_clus1_PW; %match variable name below based on this

clear y i j zz count
y=0;
count=0;
for i=1:length(xxStruct)
    if isfield(xxStruct{i}, 'm2_PW_normThetaMove')
        j=size(xxStruct{i}.m2_PW_normThetaMove,1);
        count=count+1;
        zz(count)=i;
        if length(zz)>1
            y=y+size(xxStruct{zz(count-1)}.m2_PW_normThetaMove,1);
        end
        C1_allSubj_k2_PW_normThetaWait(y+1:y+j,:)=xxStruct{i}.m2_PW_normThetaWait;
        C1_allSubj_k2_PW_normThetaMove(y+1:y+j,:)=xxStruct{i}.m2_PW_normThetaMove;
        C1_allSubj_k2_PW_normGammaWait(y+1:y+j,:)=xxStruct{i}.m2_PW_normGammaWait;
        C1_allSubj_k2_PW_normGammaMove(y+1:y+j,:)=xxStruct{i}.m2_PW_normGammaMove;
    else
        continue;
    end
end

clear y i j count zz
y=0;
count=0;
for i=1:length(xxStruct)
    if isfield(xxStruct{i}, 'nonMotor_normThetaMove')
        j=size(xxStruct{i}.nonMotor_normThetaMove,1);
        count=count+1;
        zz(count)=i;
        if length(zz)>1
            y=y+size(xxStruct{zz(count-1)}.nonMotor_normThetaMove,1);
        end
        C2_allSubj_k2_nonMotor_normThetaWait(y+1:y+j,:)=xxStruct{i}.nonMotor_normThetaWait;
        C2_allSubj_k2_nonMotor_normThetaMove(y+1:y+j,:)=xxStruct{i}.nonMotor_normThetaMove;
        C2_allSubj_k2_nonMotor_normGammaWait(y+1:y+j,:)=xxStruct{i}.nonMotor_normGammaWait;
        C2_allSubj_k2_nonMotor_normGammaMove(y+1:y+j,:)=xxStruct{i}.nonMotor_normGammaMove;
    else
        continue;
    end
end
clear y i j zz count

%% 

% compare baseline vs. exp PLV across all pairs then run ttest
sigTM=Call_allSubj_k2_motor_normThetaMove;
sigTW=Call_allSubj_k2_motor_normThetaWait;
sigGM=Call_allSubj_k2_motor_normGammaMove;
sigGW=Call_allSubj_k2_motor_normGammaWait;

count=0; %only for first subject then comment out
clear i
for i=1:size(sigTM,1)
    count=count+1;
    Call_motor_baseVec_TM(count)=mean(sigTM(i,475:1525));
    Call_motor_baseVec_TW(count)=mean(sigTW(i,475:1525));
    Call_motor_baseVec_GM(count)=mean(sigGM(i,475:1525));
    Call_motor_baseVec_GW(count)=mean(sigGW(i,475:1525));
    
    Call_motor_expVec_TM(count)=mean(sigTM(i,1550:end-500));
    Call_motor_expVec_TW(count)=mean(sigTW(i,1550:end-500));
    Call_motor_expVec_GM(count)=mean(sigGM(i,1550:end-500));
    Call_motor_expVec_GW(count)=mean(sigGW(i,1550:end-500));
end

[h,p]=ttest(baseVec_GM,baseVec_GW);


% from Drew
indB = find(t>=-2 & t<0);
indE = find(t>0 & t<=5);
for ii = 1:2
    if ii==1
        sigM = sigGM; sigW = sigGW; titl = 'gamma'; bins = -1.5:.1:1.5;
    else
        sigM = sigTM; sigW = sigTW; titl = 'theta'; bins = -3:.2:3;
    end
    WB = mean(sigW(:,indB),2); nWB = hist(WB,bins);
    WE = mean(sigW(:,indE),2); nWE = hist(WE,bins);
    MB = mean(sigM(:,indB),2); nMB = hist(MB,bins);
    ME = mean(sigM(:,indE),2); nME = hist(ME,bins);
    [~,pW] = ttest(WB,WE); [~,pM] = ttest(MB,ME);
    [~,pB] = ttest(WB,MB); [~,pE] = ttest(WE,ME);
    figure('Name',titl,'NumberTitle','off','Units','normalized','Position',[1/4 1/4 1/2 1/2],'Color','w');
    subplot(2,2,1), bar(bins,nWB,'k'); set(gca,'Box','off','XLim',[bins(1) bins(end)]); title('WAIT'); ylabel('\bf PRE'); xlabel(['p = ' num2str(pW)]);
    subplot(2,2,2), bar(bins,nMB,'k'); set(gca,'Box','off','XLim',[bins(1) bins(end)]); title('MOVE'); xlabel(['p = ' num2str(pM)]); ylabel(['p = ' num2str(pB)]);
    subplot(2,2,3), bar(bins,nWE,'k'); set(gca,'Box','off','XLim',[bins(1) bins(end)]); ylabel('\bf POST');
    subplot(2,2,4), bar(bins,nME,'k'); set(gca,'Box','off','XLim',[bins(1) bins(end)]); ylabel(['p = ' num2str(pE)]);
end



Call_allSubj_motor_normGammaWait(1:133,:)=C1_allSubj_motor_normGammaWait;
Call_allSubj_motor_normGammaWait(134:192,:)=C2_allSubj_motor_normGammaWait;
Call_allSubj_motor_normGammaWait(193:497,:)=C3_allSubj_motor_normGammaWait;

count1=0;
count2=0;
for i=1:length(clus_info)
    if clus_info(i).cluster==1
        disp('clus1')
        count1=count1+1;
        clus1chan(count1)=clus_info(i).eLbl1;
        clus1chan(count1+1)=clus_info(i).eLbl2;
        count1=count1+1;
    elseif clus_info(i).cluster==2
        disp('clus2')
        count2=count2+1;
        clus2chan(count2)=clus_info(i).eLbl1;
        clus2chan(count2+1)=clus_info(i).eLbl2;
        count2=count2+1;
    else
        disp('clus0, skipping...')
    end
end

clus1chan=unique(clus1chan);
clus2chan=unique(clus2chan);