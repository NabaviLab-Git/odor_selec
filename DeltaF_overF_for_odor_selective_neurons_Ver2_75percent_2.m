%% =======================================================================
%  Odor-selective neurons: per-cell per-trial baseline vs response
%  For each session and mouse (0L0R, 0L1R, 0L2R, 1L0R, ...)
%  Author: M. Nazari Fall 2025
% ========================================================================
clear;
close all;
clc;

%% ----- GLOBAL SETTINGS -----
mouse_tag     = '0L1R';      % <-- change to '0L1R', '0L2R', '1L0R', ...
trials_use    = 1:20;        % in case of 0L2R, thsi should be 1:18. This is the subset of trials within each session
sessions_use  = 1:8;         % session indices (Day = ses + 2)
fs            = 10;          % sampling rate (Hz)
durat         = 5;           % window length (s) for baseline & response
targetLen     = 10;          % each 10 samples of the data (1sec) is resampled to 1 sample


% What percentage of the cells should be activated during consecutive trials?
min_frac_pos = 0.75;         % 1.0 = all trials; 0.75 = most trials
min_frac_neg = 0.75;

% File naming patterns
beh_prefix    = 'MCA1';      % MCA1-<mouse_tag>_behav.csv
trace_prefix  = 'M1CA1';     % M1CA1-<mouse_tag>-Day%d-cell-traces.csv
event_prefix  = 'M1CA1';     % M1CA1-<mouse_tag> %d.csv
data_root     = mouse_tag;   % folder where traces/events live (e.g. '0L0R')

%% ----- RUN FOR NEGATIVE & POSITIVE TRIALS -----
[~, ~, results_neg] = ...
    compute_cell_activity(mouse_tag, 'neg', ...
                          beh_prefix, trace_prefix, event_prefix, data_root, ...
                          trials_use, sessions_use, fs, durat,targetLen);

[~, ~, results_pos] = ...
    compute_cell_activity(mouse_tag, 'pos', ...
                          beh_prefix, trace_prefix, event_prefix, data_root, ...
                          trials_use, sessions_use, fs, durat,targetLen);

%% ----- FIND ODOR-SELECTIVE NEURONS PER SESSION -------------------------
% For each session:
%   - pos_consistent: resp>base in ALL valid positive trials
%   - neg_consistent: resp>base in ALL valid negative trials

% per-session index lists
pos_only_idx  = cell(max(sessions_use), 1);  
neg_only_idx  = cell(max(sessions_use), 1);
both_idx      = cell(max(sessions_use), 1);
none_idx      = cell(max(sessions_use), 1);

max_ses = max(sessions_use);

% to store counts per session
n_pos_ses  = nan(max_ses,1);
n_neg_ses  = nan(max_ses,1);
n_both_ses = nan(max_ses,1);
n_none_ses = nan(max_ses,1);


for ses = sessions_use
    if ses > numel(results_pos) || ses > numel(results_neg)
        continue;
    end
    if isempty(results_pos(ses).resp_mean) || isempty(results_neg(ses).resp_mean)
        continue;
    end

    [Np, Kp] = size(results_pos(ses).resp_gt_base);
    [Nn, Kn] = size(results_neg(ses).resp_gt_base);

    if Np ~= Nn || Kp ~= Kn
        warning('Session %d: pos/neg sizes mismatch. Skipping.', ses);
        continue;
    end
    N = Np;

    pos_resp_gt = results_pos(ses).resp_gt_base;   % N x K
    pos_valid   = results_pos(ses).valid_mask;
    neg_resp_gt = results_neg(ses).resp_gt_base;   % N x K
    neg_valid   = results_neg(ses).valid_mask;

    % ---- Fraction of trials with resp > base ----
    pos_hits  = sum(pos_resp_gt & pos_valid,  2);      % # of positive trials with resp>base
    pos_total = sum(pos_valid,               2);      % # of valid positive trials
    neg_hits  = sum(neg_resp_gt & neg_valid, 2);
    neg_total = sum(neg_valid,               2);

    pos_frac = pos_hits ./ max(pos_total, 1);  % avoid division by 0
    neg_frac = neg_hits ./ max(neg_total, 1);

    pos_consistent = (pos_total > 0) & (pos_frac >= min_frac_pos);
    neg_consistent = (neg_total > 0) & (neg_frac >= min_frac_neg);

    pos_only = find( pos_consistent & ~neg_consistent );
    neg_only = find( neg_consistent & ~pos_consistent );
    both     = find( pos_consistent &  neg_consistent );
    none     = find(~pos_consistent & ~neg_consistent);

    pos_only_idx{ses} = pos_only;
    neg_only_idx{ses} = neg_only;
    both_idx{ses}     = both;
    none_idx{ses}     = none;

    % --- counts for this session 
    n_pos  = numel(pos_only);
    n_neg  = numel(neg_only);
    n_both = numel(both);
    n_none = numel(none);
    
    n_pos_ses(ses)  = n_pos;
    n_neg_ses(ses)  = n_neg;
    n_both_ses(ses) = n_both;
    n_none_ses(ses) = n_none;

    % ------ Print Info
    fprintf('\nMouse %s | Session %d (Day %d)\n', mouse_tag, ses, ses+2);
    
    fprintf('Positive-odor–selective only : %d cells (%.1f%%)\n', n_pos_ses(ses),  100 * n_pos_ses(ses)  / N);
    fprintf('Negative-odor–selective only : %d cells (%.1f%%)\n', n_neg_ses(ses),  100 * n_neg_ses(ses)  / N);
    fprintf('Both-odor responsive         : %d cells (%.1f%%)\n', n_both_ses(ses), 100 * n_both_ses(ses) / N);
    fprintf('Non-selective                : %d cells (%.1f%%)\n', n_none_ses(ses), 100 * n_none_ses(ses) / N);

   

end

% ------ Visualization
figure('Color','w');
plot(sessions_use, n_pos_ses(sessions_use),  '-o', 'LineWidth', 1.5);
xlabel('Session index');
ylabel('Number of cells');
title(sprintf('Mouse %s – cell category counts per session', mouse_tag));
grid on;
legend('Positive-selective');

figure('Color','w');
plot(sessions_use, n_neg_ses(sessions_use),  '-s', 'LineWidth', 1.5);
xlabel('Session index');
ylabel('Number of cells');
title(sprintf('Mouse %s – cell category counts per session', mouse_tag));
grid on;
legend('Negative-selective');

figure('Color','w');
plot(sessions_use, n_both_ses(sessions_use), '-^', 'LineWidth', 1.5);
xlabel('Session index');
ylabel('Number of cells');
title(sprintf('Mouse %s – cell category counts per session', mouse_tag));
grid on;
legend('Both-odor');

figure('Color','w');
plot(sessions_use, n_none_ses(sessions_use), '-d', 'LineWidth', 1.5);
xlabel('Session index');
ylabel('Number of cells');
title(sprintf('Mouse %s – cell category counts per session', mouse_tag));
legend('Non-selective');
grid on;

%% ========================================================================
%  Local function: core analysis (per mouse, per trial type)
% ========================================================================
function [sumed, percent, results] = compute_cell_activity( ...
        mouse_tag, trial_type, ...
        beh_prefix, trace_prefix, event_prefix, data_root, ...
        trials_use, sessions_use, fs, durat,targetLen)

    % trial_type: 'pos' → column 1 (positive trials)
    %             'neg' → column 2 (negative trials)

    %% ----- LOAD BEHAVIOR / TRIAL FLAGS -----
    beh_file = sprintf('%s-%s_behav.csv', beh_prefix, mouse_tag);
    if ~isfile(beh_file)
        error('Behavior file not found: %s', beh_file);
    end

    tab_beh   = readtable(beh_file);
    trl_all   = tab_beh{:,:};   % assume first two cols = [pos, neg]
    trl_pos   = trl_all(:,1);
    trl_neg   = trl_all(:,2);
    trl_pos(isnan(trl_pos)) = 0;
    trl_neg(isnan(trl_neg)) = 0;

    switch lower(trial_type)
        case 'pos'
            trial_flags = trl_pos;
        case 'neg'
            trial_flags = trl_neg;
        otherwise
            error('trial_type must be ''pos'' or ''neg''.');
    end

    num_trials_per_sess_beh = numel(trials_use);   % actual # of trials per day in behavior file
    num_trials_total        = numel(trial_flags);
    num_sessions_total      = num_trials_total / num_trials_per_sess_beh;
    if mod(num_trials_total, num_trials_per_sess_beh) ~= 0
        error('Behavior file does not have integer multiple of %d trials.', ...
              num_trials_per_sess_beh);
    end

    % trials x sessions: which trials belong to this type
    session_blocks = reshape(trial_flags, num_trials_per_sess_beh, num_sessions_total);

    %% ----- OUTPUT COLLECTORS -----
    max_ses = max(sessions_use);
    percent = cell(max_ses, 1);   % per-trial %active for each session
    sumed   = zeros(max_ses, 1);  % mean %active over selected trials
    results = struct([]);

    nRespSamples = durat * fs;
    nBaseSamples = durat * fs;
    K            = numel(trials_use);

    %% ===== MAIN SESSION LOOP =====
    for ses = sessions_use
        ses_idx    = ses;
        day_number = ses + 2;

        % Behavior: guard column existence (we index ses+1 because sessions_use starts at 1)
        if (ses+1) > size(session_blocks, 2)
            warning('Skipping ses=%d: not present in behavior blocks.', ses);
            continue;
        end

        %% ----- LOAD TRACES -----
        trace_file = fullfile(data_root, ...
            sprintf('%s-%s-Day%d-PP-BP-MC-CNMFE-cell-traces.csv', ...
                    trace_prefix, mouse_tag, day_number));
        if ~isfile(trace_file)
            warning('Trace file not found: %s. Skipping ses=%d.', trace_file, ses);
            continue;
        end

        Ttab     = readtable(trace_file);
        time_vec = Ttab{:,1}; %#ok<NASGU>
        X        = Ttab{:,2:end};      % T x N
        [T, N]   = size(X);
        X = X - mean(X);

        %% ----- LOAD EVENT TIMES -----
        event_file = fullfile(data_root, ...
            sprintf('%s-%s %d.csv', event_prefix, mouse_tag, day_number));
        if ~isfile(event_file)
            warning('Event file not found: %s. Skipping ses=%d.', event_file, ses);
            continue;
        end

        Ztab = readtable(event_file);

        % event times (sec) -> sample indices, relative to Ztab{end,2}
        z_samples_all = (Ztab(1:num_trials_per_sess_beh, 2:3).Variables - Ztab{end, 2}) * fs;
        % [trial x 2] = [onset_idx, offset_idx]

        % Alloc
        base_mean         = nan(N, K);   % mean baseline per trial
        resp_mean         = nan(N, K);   % mean response per trial
        resp_gt_base_mask = false(N, K); % response > baseline?
        valid_mask        = false(N, K); % trial used & windows valid?

        %% ----- PER-NEURON / PER-TRIAL ANALYSIS -----
        for neur = 1:N
            for k = 1:K
                trial = trials_use(k);

                % Use only trials with this trial_type (pos/neg)
                if session_blocks(trial, ses+1) == 0
                    continue;
                end

                on_off = z_samples_all(trial, :);
                if any(isnan(on_off)), continue; end

                on  = round(on_off(1));
                off = round(on_off(2));

                % Baseline window: durat s before onset
                b0 = max(1, on  - nBaseSamples + 1);
                b1 = min(T, on);

                % Response window: durat s before offset
                r0 = max(1, off - nRespSamples + 1);
                r1 = min(T, off);

                if b1 <= b0 || r1 <= r0
                    continue;
                end
                                                                                                    
                base_seg = X(b0:b1, neur);               
                resp_seg = X(r0:r1, neur);
        
                bm = mean(base_seg, 'omitnan');               
                sb = std(base_seg, 0, 'omitnan');    % baseline std (sample)


                z_base_seg = (base_seg - bm)./(sb);
                z_resp_seg = (resp_seg - bm)./(sb);

                z_base_seg_binned = resample(z_base_seg, round(durat*fs/targetLen), durat*fs);  % bining the data
                z_resp_seg_binned = resample(z_resp_seg, round(durat*fs/targetLen), durat*fs);  % bining the data

                bm = mean(z_base_seg_binned, 'omitnan');
                rm = mean(z_resp_seg_binned, 'omitnan');
                sb = std(z_base_seg_binned, 0, 'omitnan');    % baseline std (sample)


                if ~isfinite(bm) || ~isfinite(rm) || ~isfinite(sb)
                    continue;
                end
                if sb == 0
                    sb = eps;    % avoid division-by-zero / zero-threshold issues
                end
        
                base_mean(neur, k) = bm;
                resp_mean(neur, k) = rm;
                valid_mask(neur,k) = true;
        
                % CRITERION: response > baseline + 2 * baseline std, i.e.,
                % z_resp​= μ_resp​−μ_base/ std_base​​ and z_resp should be > 

                for i=1:length(z_resp_seg_binned)

                    if ((z_resp_seg_binned(i) - z_base_seg_binned(i)) / sb)> 2
                        resp_gt_base_mask(neur,k) = true;
                    end

                end
                                    
            end
        end

        fprintf('mouse %s | Day %d (ses=%d): has been analyzed\n', ...
            mouse_tag, day_number, ses);

        %% Store per-session outputs for selectivity analysis
        results(ses_idx).day          = day_number;
        results(ses_idx).base_mean    = base_mean;         % N x K
        results(ses_idx).resp_mean    = resp_mean;         % N x K
        results(ses_idx).resp_gt_base = resp_gt_base_mask; % N x K
        results(ses_idx).valid_mask   = valid_mask;        % N x K
    end
end





