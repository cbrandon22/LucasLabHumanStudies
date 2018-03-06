%% Connectivity
%
% This gets connectivity matrices from the AR1 models calculated earlier
% First threshold so that only 2 zscores stand out, then compute
% connections
%
% @author JStiso
%
%% Load Variables

top_dir = '/Users/tnl/Desktop/C/data/eeg/';
subj = 'HUP133_e';
cond = subj(end); % emergence or inductance

data_dir = [top_dir, subj, '/processed_data/dci/']; % where data is
elecs = elec_info.good;
load([data_dir, 'Analysis/order_1/AR_mod_tw300']);
save_dir = ([data_dir, 'connectivity']);
notes_dir = ['/Users/tnl/Desktop/C/data/eeg/', subj, '/docs/'];

% load anethesia start and change of consciousness
load([data_dir, 'Analysis/parameters']);
coc = numel(AR_mod) - parameters.coc; % time bin where change in consciousness occurs occurs 
ana = parameters.ana; % approximate time bin when anasthesia ws delivered/ stopped 

% calculate time period to calculate matrices over: maximum number of time
% points that will be entirely awake or entirely unconscious
dur = min([ana, coc]);

% get labels
[~, labels, ~] = xlsread([notes_dir, 'electrodeMap.xlsx']); % not sure how cell naming scheme works, so read in everything then select what you want
labels = labels(1:128,3); % should be a cell array the size of elecs when done
good_labels = labels(elecs); % should now only have labels of elecs we want to include


%% Z-score matrix

% first, create early and late connectivity matrices
A_early = zeros(size(AR_mod(1).A,1), size(AR_mod(1).A,2), dur);
A_late = zeros(size(AR_mod(1).A,1), size(AR_mod(1).A,2), dur);
B_e = zeros(size(AR_mod(1).A,1), size(AR_mod(1).A,2), dur);
B_l = zeros(size(AR_mod(1).A,1), size(AR_mod(1).A,2), dur);
S_e = zeros(size(AR_mod(1).A,1), size(AR_mod(1).A,2));
S_l = zeros(size(AR_mod(1).A,1), size(AR_mod(1).A,2));

% get data as a 3D matrix
for i = 1:dur
    A_early(:,:,i) = AR_mod(i).A;
    A_late(:,:,i) = AR_mod(end-(dur-i)).A;
end

% for each matrix, only keep values that are 2 std above the mean
for i = 1:dur
    % vectorize data
    vect_e = reshape(A_early(:,:,i),1,[]);
    vect_l = reshape(A_late(:,:,i), 1, []);
    
    % get means
    mu_e = mean(vect_e);
    mu_l = mean(vect_l);
    
    % get std
    std_e = std(vect_e);
    std_l = std(vect_l);
    
    % index
    B_e(:,:,i) = A_early(:,:,i) > mu_e + 2*std_e;
    B_l(:,:,i) = A_late(:,:,i) > mu_l + 2*std_l;
end

% get summary matrix
for i = 1:size(B_e,1)
    for j = 1:size(B_e,2)
        if i==j
            S_e(i,j) = .5;
            S_l(i,j) = .5;
        else
            S_e(i,j) = any(B_e(i,j,:));
            S_l(i,j) = any(B_l(i,j,:));
        end
    end
end

% get changes over times
S_diff = S_e - S_l;

% save
save([save_dir, '/S_e'], 'S_e', '-v7.3')
save([save_dir, '/S_l'], 'S_l', '-v7.3')
save([save_dir, '/S_diff'], 'S_diff', '-v7.3')

%% Plot

figure(1)
imagesc(S_e)
colormap('parula')
if strcmp(cond, 'i')
    title('Conscious')
else
    title('Unconscious')
end
set( gca, 'xtick', 1:numel(elecs))
set(gca, 'xticklabel', good_labels);
set(gca, 'xticklabelrotation', 60)
set( gca, 'ytick', 1:numel(elecs))
set(gca, 'yticklabel', good_labels)
colorbar

figure(2)
imagesc(S_l)
colormap('parula')
set( gca, 'xtick', 1:numel(elecs))
set(gca, 'xticklabel', good_labels);
set(gca, 'xticklabelrotation', 60)
set( gca, 'ytick', 1:numel(elecs))
set(gca, 'yticklabel', good_labels)
if strcmp(cond, 'i')
    title('Unonscious')
else
    title('Conscious')
end
colorbar

figure(3)
imagesc(S_diff)
colormap('parula')
set( gca, 'xtick', 1:numel(elecs))
set(gca, 'xticklabel', good_labels);
set(gca, 'xticklabelrotation', 60)
set( gca, 'ytick', 1:numel(elecs))
set(gca, 'yticklabel', good_labels)
title('Difference (early - late)')
colorbar

%% Better way of plotting: directed lines

% make into graph structure
g_e = digraph(S_e);
g_l = digraph(S_l);


% plot
figure
plot(g_e, 'nodelabel', good_labels)
if strcmp(cond, 'i')
    title('Conscious')
else
    title('Unconscious')
end
figure
plot(g_l, 'nodelabel', good_labels)
if strcmp(cond, 'i')
    title('Unonscious')
else
    title('Conscious')
end


