function sdat = smooth2(dat,dx,dy,varargin)
% function sdat = smooth2(dat,dx,dy,c)
%   Two-dimensional gaussian smoothing filter with standard deviations 'dx'
%   and 'dy'. A 4th input ('c') indicates that 'dat' is a circular, not
%   linear, variable (in radians). 
N = size(dat,1); T = size(dat,2);
gauw = [dx dy 5]; % gaussian smoothing kernal standard deviations and range [ms, trials, # of stds]
gwin_tim = exp(-[-gauw(1)*gauw(3)/2:gauw(1)*gauw(3)/2].^2/(2*gauw(1)^2))/(sqrt(2*pi)*gauw(1)); % gaussian for time dimension (integral = 1)
A_tim = convmtx(gwin_tim',T); lwin_tim = length(gwin_tim);
A_tim = A_tim(ceil(lwin_tim/2+1):end-floor(lwin_tim/2-1),:); % get rid of end effects
A_tim = A_tim./repmat(sum(A_tim,2),1,size(A_tim,2)); % time integral of each row is 1
gwin_trl = exp(-[-gauw(2)*gauw(3)/2:gauw(2)*gauw(3)/2].^2/(2*gauw(2)^2))/(sqrt(2*pi)*gauw(2)); % gaussian for trial dimension (integral = 1)
A_trl = convmtx(gwin_trl',N); lwin_trl = length(gwin_trl);
A_trl = A_trl(ceil(lwin_trl/2+1):end-floor(lwin_trl/2-1),:); % get rid of end effects
A_trl = A_trl./repmat(sum(A_trl,2),1,size(A_trl,2)); % time integral of each row is 1
if nargin < 4 % 'dat' is linear variable %%% TODO: get rid of repmat %%%
    if any(isnan(dat(:))), for ist = 1:T, sdat(:,ist) = nansum(A_trl.*repmat(dat(:,ist)',N,1),2); end % slow, but robust to NaNs
    else, for ist = 1:T, sdat(:,ist) = A_trl*dat(:,ist); end, end
    for ist = 1:N, sdat(ist,:) = (A_tim*sdat(ist,:)')'; end
else % 'dat' is circular variable (radians) %%% TODO: add NaN support %%%
    for ist = 1:T, C = A_trl*cos(dat(:,ist));       S = A_trl*sin(dat(:,ist));      sdat(:,ist) = atan2(S,C); end
    for ist = 1:N, C = (A_tim*cos(sdat(ist,:)'))';  S = (A_tim*sin(sdat(ist,:)'))'; sdat(ist,:) = atan2(S,C); end
end