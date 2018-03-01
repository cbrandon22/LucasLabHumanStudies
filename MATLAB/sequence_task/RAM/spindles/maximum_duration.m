function [bv,begins,ends] = maximum_duration(bv,begins,ends,max_dur,fs)
        % MAXIMUM_DURATION - checks the sample duration of the spindles/ripples.
        % Input is a vector containing ones in the interval where the spindle/ripple is
        % and indexs describing the start and end of the spindle/ripple. The last two
        % inputs are the maximum duration given in seconds and the sampling
        % frequency given in Hz.
        % Output is a vector containing ones in the interval where the spindle/ripple with
        % duration shorter than or equal to the maximum duration is and indexs
        % describing the start and end of the spindle/ripple.
        
        duration_samples = ends-begins+1;
        for k = 1:length(begins)
            if duration_samples(k) > max_dur*fs
                bv(begins(k):ends(k)) = 0;
                begins(k) = 0;
                ends(k) = 0;
            end
        end
        begins = begins(begins~=0);
        ends = ends(ends~=0);
    end