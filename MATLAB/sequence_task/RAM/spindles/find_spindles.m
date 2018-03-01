function [begins, ends] = find_spindles(bv)
        % FIND_SPINDLES - find start and end index' of spindles.
        % Input is a binary vector bv containing ones where spindles (or ripples) are detected.
        % Output is vectors containing the index' of spindle/ripple beginnings and ends
        % (first sample of spindle/ripple and last sample of spindle/ripple, respectively).
        
        sise = size(bv);
        E = bv(2:end)-bv(1:end-1); % Find start and end of intervals with spindles
        
        begins = find(E==1)+1;
        if bv(1) == 1
            if sise(1) > 1
                begins = [1; begins];
            elseif sise(2) > 1
                begins = [1 begins];
            else
                error('The input signal is not one dimensional')
            end
        elseif numel(begins) == 0 && bv(1) == 0
            begins = NaN;
        end
        
        ends = find(E==-1);
        if bv(end) == 1
            if sise(1) > 1
                ends = [ends; length(bv)];
            elseif sise(2) > 1
                ends = [ends length(bv)];
            else
                error('The input signal is not one dimensional')
            end
        elseif numel(ends) == 0 && bv(end) == 0
            ends = NaN;
        end
    end

    

    