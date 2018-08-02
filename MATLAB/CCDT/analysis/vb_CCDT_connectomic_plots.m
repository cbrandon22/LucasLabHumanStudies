%%
% load Wst0 (W_st vRT vDT chlbl)
ploton = 1;

%%
Wprecue = zeros(length(vRT),length(chlbl),length(chlbl),3);
WDTone = zeros(length(vRT),length(chlbl),length(chlbl),3);
WDTtwo = zeros(length(vRT),length(chlbl),length(chlbl),3);
WDTthree = zeros(length(vRT),length(chlbl),length(chlbl),3);

%%
for freq=1:3
    for i=1:length(vRT)
        Wprecue(i,:,:,freq) = W_st(i).precue(freq).adj;
        WDTone(i,:,:,freq) = W_st(i).DTone(freq).adj;
        WDTtwo(i,:,:,freq) = W_st(i).DTtwo(freq).adj;
        WDTthree(i,:,:,freq) = W_st(i).DTthree(freq).adj;
    end
end

%% extract from W
nR_precue = zeros(length(vRT),length(chlbl),3);
nR_DTone = zeros(length(vRT),length(chlbl),3);
nR_DTtwo = zeros(length(vRT),length(chlbl),3);
nR_DTthree = zeros(length(vRT),length(chlbl),3);

cC_precue = zeros(length(vRT),length(chlbl),3);
cC_DTone = zeros(length(vRT),length(chlbl),3);
cC_DTtwo = zeros(length(vRT),length(chlbl),3);
cC_DTthree = zeros(length(vRT),length(chlbl),3);



for i=1:length(vRT)
    for fq = 1:3
        str_precue(i,:,fq) = strengths_und(squeeze(Wprecue(i,:,:,fq)));
        %         nR_precue(i,:,fq) = betweenness_wei(squeeze(Wprecue(i,:,:,fq)));
        netQexp_precue(:,:,i,fq) = getCommunicability(squeeze(Wprecue(i,:,:,fq)),1,1);
        %         cC(i,:,fq) = clustering_coef_wu(squeeze(Wprecue(i,:,:,fq)));
        %         [~,D] = reachdist(squeeze(Wprecue(i,:,:,fq)));
        %         pL(i,:,fq) = mean(D,1);
        %         clear D
        %         bC
    end
end
nodalQexp = squeeze(mean(netQexp_precue,1));
nodalQexp = permute(nodalQexp,[2 1 3]);
netQexp_precue = permute(netQexp_precue,[3 1 2 4]);

%%
%finding significant predictors
% assess each metric for trial-by-trial predictive correlation
% for freq=1:3
% for i=1:length(chlbl)
% [r_cMRT(i), p_cMRT(i)] = corr(nodalQexp(vRT>0&vDT>1450,i,freq),vRT(vRT>0&vDT>1450));
% end
% plot the trend for verification
% figure;
% scatter(r_cMRT, p_cMRT)
% % sort by p-values to find most predictive
% [pV,pI] = sort(p_cMRT, 'ascend');
% pV(1:10)
% pI(1:10) %channels
% % sort by correlations
% [rV,rI] = sort(r_cMRT ,'ascend');
% %check
% pI(1)-rI(1); %should be the same, thus should equal zero

% %rename to store
% r_Qexp(freq) = r_cMRT;
% p_Qexp(freq) = p_cMRT;

% end

%%
for ix = 1:3
    switch ix
        case 1
            currDT = 1500;
            currTrials = (vRT>0&vDT>1450);
            
            [a,b] = sort(vRT, 'ascend');
            cTc = currTrials;
            cB = (b(cTc(b)));
            iF = zeros(length(vRT),1);
            iS = zeros(length(vRT),1);
            iE = zeros(length(vRT),1);
            iL = zeros(length(vRT),1);
            iF(cB(1:round(sum(cTc)/3))) = 1;
            iS(cB(end-round(sum(cTc)/3)+1:end)) = 1;
            iE(1:round(sum(cTc)/3)) = 1;
            iL(end-round(sum(cTc)/3)+1:end) = 1;
            currFi = logical(iF);
            currSi = logical(iS);
            currEi = logical(iE);
            currLi = logical(iL);
            sumChk1500 = [sum(currFi) sum(currSi) sum(currLi) sum(currEi)];
            if diff(sumChk1500)==[0 0 0]
               ifast1500 = currFi;
               islow1500 = currSi;
               ilate1500 = currLi;
               iearly1500 = currEi;
            else
                disp('error with group indices')
                break
            end

            
        case 2
            currDT = 500;
            currTrials = (vRT>0&vDT>450&vDT<550);
            
            [a,b] = sort(vRT, 'ascend');
            cTc = currTrials;
            cB = (b(cTc(b)));
            iF = zeros(length(vRT),1);
            iS = zeros(length(vRT),1);
            iE = zeros(length(vRT),1);
            iL = zeros(length(vRT),1);
            iF(cB(1:round(sum(cTc)/3))) = 1;
            iS(cB(end-round(sum(cTc)/3)+1:end)) = 1;
            iE(1:round(sum(cTc)/3)) = 1;
            iL(end-round(sum(cTc)/3)+1:end) = 1;
            currFi = logical(iF);
            currSi = logical(iS);
            currEi = logical(iE);
            currLi = logical(iL);
            sumChk500 = [sum(currFi) sum(currSi) sum(currLi) sum(currEi)];
            if diff(sumChk500)==[0 0 0]
               ifast500 = currFi;
               islow500 = currSi;
               ilate500 = currLi;
               iearly500 = currEi;
            else
                disp('error with group indices')
                break
            end
            
        case 3
            currDT = 0;
            guesses = (vRT<200&vRT>-200);
            errors = (vRT<-200);
            currTrials = (vRT<200);
            
    end
    
    for x=1:3 %1=theta, 2=high gamma, 3=low gamma
        if ploton
            figure;
            subplot(2,3,1)
            plot(vRT, 'k.')
            hold on
            if ix==1||ix==2
                plot(find(currSi),vRT(currSi), 'mo')
                plot(find(currFi),vRT(currFi), 'go')
                title(['Behavior (RT) DT: ' num2str(currDT) 'ms'] )
                ylim([0 max(vRT)])
            else
                plot(find(guesses),vRT(guesses),'go')
                plot(find(errors),vRT(errors),'mo')
                title('Behavior (RT) 0+/-200ms vs. Errors')
                ylim([min(vRT) 200])
            end
            subplot(2,3,2)
            if ix==1||ix==2
                plot(mean(str_precue(currFi,:,x),1), 'g-')
                hold on
                plot(mean(str_precue(currSi,:,x),1), 'm-')
            else
                plot(mean(str_precue(guesses,:,x),1),'g-')
                hold on
                plot(mean(str_precue(errors,:,x),1),'m-')
            end
            box off
            if x==1
                title('3-12Hz FS nodal strength')
            elseif x==2
                title('70-100Hz FS nodal strength')
            else
                title('35-55Hz FS nodal strength')
            end
            subplot(2,3,3)
            if ix==1||ix==2
                plot(mean(str_precue(currLi,:,x),1), 'b-')
                hold on
                plot(mean(str_precue(currEi,:,x),1), 'k-')
                box off
                if x==1
                    title('3-12Hz LE nodal strength')
                elseif x==2
                    title('70-100Hz LE nodal strength')
                else
                    title('35-55Hz LE nodal strength')
                end
            elseif ix==3&&sum(guesses)>0
                scatter(mean(str_precue(guesses,:,x),2), vRT(guesses), 'ko')
                currCM = mean(str_precue(guesses,:,x),2);
                pcMRT=polyfit(currCM,vRT(guesses),1);
                fitcMRT=polyval(pcMRT,currCM);
                [r_cMRT, p_cMRT] = corr(currCM,vRT(guesses));
                hold on
                plot(currCM,fitcMRT,'r--')
                text(max(get(gca,'XLim')),max(get(gca,'YLim')),['r = ' num2str(r_cMRT,'%4.3f') ''],'HorizontalAlignment','right','VerticalAlignment','top', 'color', 'black');
                text(max(get(gca,'XLim')),max(get(gca,'YLim')),['p = ' num2str(p_cMRT,'%4.2f')],'HorizontalAlignment','left','VerticalAlignment','top', 'color', 'red');
                box off
                if x==1
                    title('0+/-200 only RTf(3-12Hz netStrength)')
                elseif x==2
                    title('0+/-200 only RTf(70-100Hz netStrength)')
                else
                    title('0+/-200 only RTf(35-55Hz netStrength)')
                end
                set(gcf,'color','w')
                
            end
            
            subplot(2,3,4)
            scatter(mean(str_precue(currTrials,:,x),2), vRT(currTrials), 'ko')
            currCM = mean(str_precue(currTrials,:,x),2);
            pcMRT=polyfit(currCM,vRT(currTrials),1);
            fitcMRT=polyval(pcMRT,currCM);
            [r_cMRT, p_cMRT] = corr(currCM,vRT(currTrials));
            hold on
            plot(currCM,fitcMRT,'r--')
            text(max(get(gca,'XLim')),max(get(gca,'YLim')),['r = ' num2str(r_cMRT,'%4.3f') ''],'HorizontalAlignment','right','VerticalAlignment','top', 'color', 'black');
            text(max(get(gca,'XLim')),max(get(gca,'YLim')),['p = ' num2str(p_cMRT,'%4.2f')],'HorizontalAlignment','left','VerticalAlignment','top', 'color', 'red');
            box off
            if x==1
                title('RTf(3-12Hz netStrength)')
            elseif x==2
                title('RTf(70-100Hz netStrength)')
            else
                title('RTf(35-55Hz netStrength)')
            end
            set(gcf,'color','w')
            
            subplot(2,3,5)
            allCM = nodalQexp(currTrials,:,x);
            [r_cMRT_all, p_cMRT_all] = corr(allCM,vRT(currTrials));
            [~,maxSig] = min(p_cMRT_all);
            
            scatter(nodalQexp(currTrials,maxSig,x),vRT(currTrials))
            hold on
            currCM = nodalQexp(currTrials,maxSig,x);
            pcMRT=polyfit(currCM,vRT(currTrials),1);
            fitcMRT=polyval(pcMRT,currCM);
            [r_cMRT, p_cMRT] = corr(currCM,vRT(currTrials));
            hold on
            plot(currCM,fitcMRT,'r--')
            text(max(get(gca,'XLim')),max(get(gca,'YLim')),['r = ' num2str(r_cMRT,'%4.3f') ''],'HorizontalAlignment','right','VerticalAlignment','top', 'color', 'black');
            text(max(get(gca,'XLim')),max(get(gca,'YLim')),['p = ' num2str(p_cMRT,'%4.2f')],'HorizontalAlignment','left','VerticalAlignment','top', 'color', 'red');
            box off
            if x==1
                title(['3-12Hz Single-Trial Qexp: Optimal Ch ' num2str(maxSig)])
            elseif x==2
                title(['70-100Hz Single-Trial Qexp: Optimal Ch ' num2str(maxSig)])
            else
                title(['35-55Hz Single-Trial Qexp: Optimal Ch ' num2str(maxSig)])
            end
            set(gcf,'color','w')
            
        end
    end
end


