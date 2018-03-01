function create_surr_data_bl(dir, surr_iter, win_start, win_end, elecs, NumEvents,bl_start,bl_end,ictal_events)
%
%   create_surr_data(dir, surr_iter, win_start, win_end, elecs, NumEvents [, bl_start, bl_end, ictal_events])


    if (nargin < 9)
        ictal_events = [];
    end
    elecs  = sort(elecs);
    tm = 0; 
    total = tic;
    for e = elecs
        load([dir 'rawSPEC_e' int2str(e)]);
        dat_srate = rawSPEC.srate;
        dat_length = length(rawSPEC.SPEC).*rawSPEC.compression;
        win_length = ceil(win_end/1000*dat_srate/rawSPEC.compression)-floor(win_start/1000*dat_srate/rawSPEC.compression)+1;
        bl_length = ceil(bl_end/1000*dat_srate/rawSPEC.compression)-floor(bl_start/1000*dat_srate/rawSPEC.compression)+1;
         % if this the first iteration create random correct onsets that do
         % not fall on ictal periods as defined by ictal events
        if (find(elecs==e)==1)   
            onsets = 1:dat_length;
            for i = 1:size(ictal_events,1)
                start_ictal = floor(ictal_events(i,1)*dat_srate)-(ceil(win_end/1000*dat_srate)-floor(bl_start/1000*dat_srate)+1);    % start of invalid win_start point
                end_ictal   = min(ceil(ictal_events(i,2)*dat_srate),dat_length);                                                      % end of invalid win_end point
                if (start_ictal < 1)
                    start_ictal = 1;
                end
                onsets((onsets>= start_ictal & onsets<= end_ictal)) = [];
            end
            
            if (win_start<0 || bl_start<0)
                fl_start = min(win_start,bl_start);
                onsets((onsets<= (floor(abs(fl_start)/1000*dat_srate)+1+rawSPEC.compression))) = [];  % remove the first time stamps (not to go before the first time point)
            end
            onsets((onsets>= (dat_length-ceil(win_end/1000*dat_srate)-1-rawSPEC.compression))) = [];% remove the last time stamps (not to go beyond the last time point)
        end

         % Print out time estimate
        est = (length(elecs)-find(elecs==e)+3)*(tm/(find(elecs==e)-1));
        fprintf('Electrode %d of %d. Estimated time left %d:%02d\n',e, elecs(end),floor(est/60),round(rem(est,60)));
        tic;

        surr_distrib = zeros(length(rawSPEC.freqs),win_length);
        surr_distrib_sq = zeros(length(rawSPEC.freqs),win_length);
        for i=1:surr_iter
            surr_onsets = onsets(randi(length(onsets),[1 NumEvents]));          % pool NumEvents random time stamps from onsets array
            try
                tmp = create_spec_erp_fast_bl(rawSPEC,surr_onsets,win_start,win_end,surr_onsets,bl_start,bl_end);
            catch ME
                disp('Error, please check !!!!');
                keyboard
            end
            surr_distrib(:,:)= surr_distrib(:,:) + tmp;
            surr_distrib_sq(:,:) = surr_distrib_sq(:,:) + tmp.^2;
        end
        tm = tm+toc;
        MN(:,:,e) = surr_distrib./surr_iter;
        SD(:,:,e) = sqrt(...
                         (surr_distrib_sq - (surr_distrib).^2./surr_iter)... 
                        ./(surr_iter-1));
    end
    freqs = rawSPEC.freqs;
    save([dir 'erps_surr_data_bl_',num2str(elecs(1)),'_' num2str(NumEvents) 'evs_' num2str(win_start) '_' num2str(win_end) '_' num2str(bl_end)],'MN','SD','win_start','win_end','freqs','bl_start','bl_end');
    fprintf('Total elapsed time: %4.3f minutes\n',toc(total)./60);
end
