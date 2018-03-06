function [ min_array ] = get_min_interval( vect )
%Gets the minimum interval of the input array- useful for finding bin edges

% initial guess at min
curr_min = inf;
sort_array = sort(vect);
for i = 2:numel(sort_array)
    min_guess = abs(sort_array(i-1)- sort_array(i));
    if curr_min > min_guess && min_guess ~= 0;
        curr_min = min_guess;
    end
end
min_array = curr_min;

end

