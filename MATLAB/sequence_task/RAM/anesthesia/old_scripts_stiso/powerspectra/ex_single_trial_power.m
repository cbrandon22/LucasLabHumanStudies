%% define variables to run
% This script will plot 64 electrodes per figure. If you want to plot only

subj = 'ST43'; % Define for saving purposes, sometimes should be with hand e.g. GP48/Left
condition = 'GraspL';

f1 = 70; %lower bound of gamma
f2 = 150; %upper bound of gamma

start_time_window = -400; % values in milliseconds, -1200 for view, -400 to grasp, 0 for rest (change in BL below too)
end_time_window = 5000; %2000 for view, 5000 for grasp 2000 for reset 
%change according to # of elecs

e1=[1:64];
e_trans1=[1:64] ;

factor = 0; % Changed if plotting more than 64 electrodes. 
            % ie: for electrodes 1-64 factor = 0, 65-128, factor = 1, 129-192 factor = 2, etc. 
            % If plotting fewer than 64 electrodes, factor only effects what the name of the picture is. 
            % In other words, factor is only used to define plot elecs below, and in saving.

% set dimension of divide screen. Works best if dimensions are square.
% remember to change name of jpeg created (line 118). If plotting only one
% elec change to 1X1
width = 8;
height = 8;

plot_elecs = [64:-1:1]+(64*factor); % electrodes to plot; 
                                 % only change if also changing dimensions of divide screen. 
                                 % For example, if you want to only plot one electrode, set plot_elecs to the index 
                                 % of that electrode, and both dimensions to 1

%% loading variables

load(['/Users/Anat/Documents/Berkeley lab/MNS hands ECoG/',subj,'/srate'])

load(['/Users/Anat/Documents/Berkeley lab/MNS hands ECoG/',subj,'/gdat_CAR_good'])

% load(['/Users/Anat/Documents/Berkeley lab/MNS hands ECoG/',subj,'/onsetsGraspR'])
% load(['/Users/Anat/Documents/Berkeley lab/MNS hands ECoG/',subj,'/onsetsRestR'])
% load(['/Users/Anat/Documents/Berkeley lab/MNS hands ECoG/',subj,'/onsetsViewR'])

load(['/Users/Anat/Documents/Berkeley lab/MNS hands ECoG/',subj,'/onsetsGraspL'])
load(['/Users/Anat/Documents/Berkeley lab/MNS hands ECoG/',subj,'/onsetsRestL'])
load(['/Users/Anat/Documents/Berkeley lab/MNS hands ECoG/',subj,'/onsetsViewL'])

eval (['onsets=onsets' condition ';'])  


%setting up windows and timing


plot_jump = 200;
tm_st = round(start_time_window./1000*srate);
tm_en = round(end_time_window./1000*srate);
jm = round(plot_jump./1000*srate);

bl_st = round(-5200 ... %change BL accordingly. Now changed to BL of view for all conditions.400 ms starting at -600 for View, -5200 for grasp and -2800 for rest  
    ./1000*srate);
bl_en = round(-4800  ...
    ./1000*srate);


clr_res = 'k';


%setting up the figure (need DivideScreen script)

ScSz=[1 1 1680 1050];
pos_grid=[1 1 1680 1050];
fgrid = figure;
set(fgrid, 'Position', ScSz,'color', [1 1 1], 'MenuBar', 'figure');
gridPos = DivideScreen(height, width, ScSz, 30, 60, pos_grid);


% Make graphs


e_titles = {'grid'};

%plotting each electrode

for elec = plot_elecs %Change according to what you want to plot
    
    %extracting gamma analytic amplitude for each electrode
    
    band = gdat_CAR_good(elec, :);
    band = eegfilt(band,srate,[],f2);
    band = eegfilt(band,srate,f1,[]);
    band = abs(hilbert(band));
    band =  band_pass(band, srate, 0.1, 10); 
    name='Gamma';
    
    clear TrialsMTX; clear tm_stmps; clear bl_stmps;
    
    %setting up onsets for each trial and creating data matrix of trials x
    %time points, baseline corrected and in %change of baseline
    
    for i = 1:(length(onsets))
        tm_stmps = (onsets(i)+tm_st) : (onsets(i)+tm_en);
        bl_stmps = (onsets(i)+bl_st) : (onsets(i)+bl_en);
        TrialsMTX(i,:) = 100 * ((band(tm_stmps) - mean(band(bl_stmps), 2)))./mean(band(bl_stmps), 2); % percent change from baseline
    end
    
    if ~exist(['/Users/Anat/Documents/Berkeley lab/MNS hands ECoG/',subj, '/Analysis/',condition,], 'dir')
        mkdir(['/Users/Anat/Documents/Berkeley lab/MNS hands ECoG/',subj], ['/Analysis/', condition]);
    end
    
    if ~exist(['/Users/Anat/Documents/Berkeley lab/MNS hands ECoG/',subj,'/Analysis/',condition,'singleTrials_HG_400viewBL'], 'dir') %all conditions base lined to -400->-200 before view
        mkdir(['/Users/Anat/Documents/Berkeley lab/MNS hands ECoG/',subj, '/Analysis/', condition,], ['singleTrials_HG_400viewBL']);
    end
    save(['/Users/Anat/Documents/Berkeley lab/MNS hands ECoG/',subj,'/Analysis/',condition,'/singleTrials_HG_400viewBL/TrialsMTX_e' num2str(elec)],'TrialsMTX')    
    
   % set axis labels based on condition
   switch upper(condition(1))
       case 'G'
           ticks = [400,2400,4400];
           labels = {'0','2000','4000'};
       case 'V'
           ticks = [200,1200,2200];
           labels = {'-1000','0','1000'};
       case 'R'
           ticks = [0,500,1000,1500];
           labels = {'0', '500', '1000','1500'};
   end
    
    %actual plotting of each electrode
    h=axes('Position',gridPos{e_trans1(find(elec == plot_elecs))})%elec-(64*(factor)))}); %factor changes for each group of elecs we plot
    
    
    pcolor(double(TrialsMTX));
    %pcolor(tm_st:(tm_en-1),1:i,double(TrialsMTX(:,:)'));
    shading flat;
    caxis([-150 150]);
    title(sprintf('%s %d',name,elec));
    %xlabel('ms'); ylabel('Trial #');
    for z = 1:length([tm_st:jm:tm_en])
        plot_str{z} = start_time_window+(z-1)*plot_jump;
    end
    hold on;
    set(gca, 'XTick',ticks);
    set(gca, 'XTickLabel', labels);
    plot([abs(start_time_window),abs(start_time_window)], [1 numel(onsets)],'k', 'LineWidth', 2); %should plot black line at 0 (parallel to y axis)
   
    clear band; 
    set(fgrid,'PaperPositionMode','auto')
    
end
colorbar;
saveas(gcf, ['/Users/Anat/Documents/Berkeley lab/MNS hands ECoG/',subj,'/Analysis/',condition,'/singleTrials_HG_400viewBL/single_trial_power_',num2str(factor+1),'.jpg'], 'jpg') %changes folder name occurding to what you're saving

%close(fgrid)

%%
% 
% bl_st = round( -200 ./1000*srate);
% bl_en = round( 0 ./1000*srate);
% start_time_window = -200;
% tm_st = round( start_time_window ./1000*srate);
% tm_en = round( 5000 ./1000*srate);
% plot_jump = 200;
% jm = round(plot_jump./1000*srate);
% f1=70;
% f2=150;
% lines(1,1:5201)=0; %should be the same "length" as TrialsMTX_eN - to create the red line at zero
% % 
% % 
% % for i=1:64
% %     
% %     load(['TrialsMTX_e' num2str(i)])
%     
%     fgrid=figure;
%     
%     plot(tm_st:tm_en, mean(TrialsMTX), 'LineWidth', 3);
%   
%     title(sprintf('%d, Power %d->%d %change',elec, f1, f2));
%         
%     xlabel('ms'); 
%     ylabel('%change');
%     
%     hold on;
%     plot(tm_st:tm_en, lines, 'r', 'LineWidth', 3);
% 	
%     
%     for z = 1:length([tm_st:jm:tm_en])
%         plot_str{z} = start_time_window+(z-1)*plot_jump; 
%     end
%     
%     set(gca, 'XTick', [tm_st:jm:tm_en], 'XTickLabel', plot_str, 'XTickMode', 'manual', 'Layer', 'top');
%     xlim([tm_st tm_en])
%     
%     hold off;
%     
%     set(fgrid, 'PaperPositionMode', 'auto')
%    
%     print(fgrid, '-djpeg', ['/Users/Anat/Documents/Berkeley lab/MNS hands ECoG/GP48/Left/Analysis/GraspL/singleTrials_HG/Mean_trace_e' num2str(i) '.jpg'])
%  
%    close(fgrid);
%     
% end

% %% Stats
% 
% 
% for i=1:64
% load(['TrialsMTX_e' num2str(i)])
% [h, p, ci, stats]=ttest(mean(TrialsMTX(:,round(srate):round((2.250*srate))),2)); %750 to 2000 ms
% active_stats.motor.pvalue(i,1)=p;
% active_stats.motor.tvalue(i,1)=stats.tstat;
% end
% 
% for i=1:64
% load(['TrialsMTX_e' num2str(i)])
% [h, p, ci, stats]=ttest(mean(TrialsMTX(:,(round(.5*srate)):round(srate)),2)); %250 ms to 750 ms - plus the baseline, which is in the TrialsMTX file
% active_stats.cue.pvalue(i,1)=p;
% active_stats.cue.tvalue(i,1)=stats.tstat;
% end
% 
% save('active_stats','active_stats')
% 
% 
% 
