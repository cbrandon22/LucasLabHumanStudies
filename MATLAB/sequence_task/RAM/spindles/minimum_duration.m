function [bv,begins,ends] = minimum_duration(bv,begins,ends,min_dur,fs)
        % MINIMUM_DURATION - checks the sample duration of the spindles or ripples.
        % Input is a vector containing ones in the interval where the spindle/ripple is
        % and indexs describing the start and end of the spindle/ripple. The last two
        % inputs are the minimum duration given in seconds and the sampling
        % frequency given in Hz.
        % Output is a vector containing ones in the interval where the spindle/ripple with
        % duration longer than or equal to the minimum duration is and indexs
        % describing the start and end of the spindle.
        
        duration_samples = ends-begins+1;
        for k = 1:length(begins)
            if duration_samples(k) < min_dur*fs
                bv(begins(k):ends(k)) = 0;
                begins(k) = 0;
                ends(k) = 0;
            end
        end
        begins = begins(begins~=0);
        ends = ends(ends~=0);
    end