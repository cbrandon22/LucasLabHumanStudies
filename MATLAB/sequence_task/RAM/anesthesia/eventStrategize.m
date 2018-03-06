function [correctR incorrectR] = eventStrategize(events)
ci = 1;
ii = 1;
% sr = 2.047999e+03 %sample rate from 


for i = 1:length(events)
    if strcmp(events(i).type,'RECALL') && strcmp(events(i).response,events(i).target);
        correctR(ci) = events(i);
        ci = ci + 1;
    elseif strcmp(events(i).type,'RECALL')
        incorrectR(ii) = events(i);
        ii = ii + 1;

    end
end