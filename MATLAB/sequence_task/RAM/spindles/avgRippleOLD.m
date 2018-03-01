function avgSpindle
%         win1 = 500;
%         win2 = 500;
%         figure('Name','EEG Spindle')
        i = 1;
    for j = 1:length(detected_spindles)
        probs(j) = detected_spindles{1,j}.prob_spindle;
        duration(j) = detected_spindles{1,j}.spindle_duration_sec;
        Spind = j;
        
        if detected_spindles{1,j}.prob_spindle > 0.4; % spindle probability threshold
            
            SP(i) = Spind;
%             V = dat(spindles_start_end(Spind,1)-win1:spindles_start_end(Spind,2)+win2);
%             T = t(spindles_start_end(Spind,1)-win1:spindles_start_end(Spind,2)+win2);
%             T2 = 0-win2:spindles_start_end(Spind,2)+win2-spindles_start_end(Spind,1);
%             VZ = zscore(V);
%             Spindle{1,i}.VZ = VZ;
%             Spindle{1,i}.T2 = T2;
%             SpindleVZ{i} = VZ;
%             SpindleT2{i} = T2;
%            plot(T2/1e3,VZ);
%            hold on
           
           i = i+1;
        end
    end

%     for ich = 1:size(PM,1)
        PM(ich,:,:) = smooth2(squeeze(PM(ich,:,:)),2,1);
        PM(ich,:,:) = zscore(squeeze(PM(ich,:,:))')'; % zscore
%     end

% plot Spectrogram
        ich = 1;
        tank = 'ditka_cage_040615'; % tank 
%         for ich = 1:size(PM,1)
            figure('Name',[tank ': ' num2str(ich)],'NumberTitle','off','Units','normalized','Position',[1/8 1/4 3/4 1/2],'Color','w');
            imagesc(x,1:length(freq),squeeze(PM(ich,:,:)));
            ytick = [1:10 20:10:100 200];
            yticknm = {'1';'';'';'';'';'';'';'';'';'10';'';'';'';'';'';'';'';'';'100';''};
            yind = zeros(size(ytick));
            for ii = 1:length(ytick)
                [~,yind(ii)] = min(abs(freq-ytick(ii)));
            end
            set(gca,'Box','off','YDir','normal','TickDir','out','YLim',[1 length(freq)],'XLim',[x(1) x(end)],'YTick',yind,'YTickLabel',yticknm);
            xlabel('time (h)'); ylabel('frequency (Hz)'); title([tank ': ' num2str(ich)],'Interpreter','none');
            set(gca,'CLim',[-3 5]);
            hc = colorbar; set(get(hc,'YLabel'),'String','power');
%         end
        
      
%         oness = ones(1,length(spindles_start_end(SP,1)));
%         oness = oness*10;
%         PMa = mean(PMs);

        % Calculate Instantaneous Frequencies for the entire dataset
        fs = 1000;
        datastart = 1;
        if datastart == 1
            timecorrection = 1; % convert to allow hours entry above, use 1 if datastart = 1;
        else 
            timecorrection = 3600*fs; % convert to allow hours entry above, use 1 if datastart = 1;
        end
        dataend = 10;
        data = dat(datastart*timecorrection:dataend*3600*fs);
        [b,a] = butter(2,[0.5 16]/(fs/2)); % delta & theta band filter 0.5-16
        data = filtfilt(b,a,data);
%         delta_amp = abs(hilbert(data)); % delta band amplitude
%         z = hilbert();
%         instamp = abs(z);
        instfreq = fs/(2*pi)*diff(unwrap(angle(hilbert(data)))); %Compute the analytic signal and differentiate its phase to measure the instantaneous frequency. The scaled derivative yields a meaningful estimate.
    
%         xTime = find(x==datastart);
%         xTime2 = find(x==dataend);
%         PMs = squeeze(PM(:,:,:));
%         xFreq = PMs(:,xTime:xTime2);
%         xTime3 = x(xTime:xTime2);
%         plot(xTime3,xFreq)
%         hold on
        
%         surf(x(xTime:xTime2),freq,PMs(:,xTime:xTime2))
%         hold on
        
        % Take Instantaneous frequencies of each spindle start time
        if timecorrection == 1
            datastart1 = 0;
        else
            datastart1 = datastart;
        end
        iFreq = instfreq(spindles_start_end(SP,1));
        spindT = spindles_start_end(SP,1)/3600/1000.+datastart1; % spindle start times in dimensions of x (hours), shifted to the start of spindle data
        hold on
        scatter(spindT,abs(iFreq))
        hold off
        figure('Name','Spindle Count')
        histogram(spindles_start_end(:,1)/3600/1000.+datastart1, 'DisplayStyle','stairs')
        % Create xlabel
        xlabel({'Time (h)'});

        % Create ylabel
        ylabel({'Spindle Count'});

        % Create title
        title({'Spindle Counts'});
        
%         freqP = PMs(:,spindT);
%         stairs(spindT,PM(spindT))
%         bar(spindT,oness)
        
    % K cluster
    
%     [IDX, C, SUMD] = kmeans(squeeze(PM), 4);
    
    [IDX3, C3, SUMD3] = kmeans([iFreq spindT], 4);
    STFK = [spindT iFreq IDX3];
    STFKsort = sortrows(STFK,3);
%     scatter(STFKsort(:,1),STFKsort(:,3).*5) % plot time x kmean cluster
    
    scatter(STFKsort(STFKsort(:,3)==1,1),STFKsort(STFKsort(:,3)==1,2))
    hold on
    scatter(STFKsort(STFKsort(:,3)==2,1),STFKsort(STFKsort(:,3)==2,2))
    hold on
    scatter(STFKsort(STFKsort(:,3)==3,1),STFKsort(STFKsort(:,3)==3,2))
    hold on
    scatter(STFKsort(STFKsort(:,3)==4,1),STFKsort(STFKsort(:,3)==4,2))
    
    % Plot on spectrogram - ruins the grouping.
    
    hold on    
    scatter(STFKsort(STFKsort(:,3)==1,1),abs(STFKsort(STFKsort(:,3)==1,2))*10)
    hold on
    scatter(STFKsort(STFKsort(:,3)==2,1),abs(STFKsort(STFKsort(:,3)==2,2))*10)
    hold on
    scatter(STFKsort(STFKsort(:,3)==3,1),abs(STFKsort(STFKsort(:,3)==3,2))*10)
    hold on
    scatter(STFKsort(STFKsort(:,3)==4,1),abs(STFKsort(STFKsort(:,3)==4,2))*10)
    
end