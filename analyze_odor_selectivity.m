%% ========================================================================
%  ODOR SELECTIVITY ANALYSIS
%  
%  Identifies odor-selective neurons from calcium imaging data across
%  behavioral training sessions. Classifies neurons as selective to
%  positive odors, negative odors, both, or non-selective.
%
%  Author: M. Nazari
%  Date: Fall 2025
%  
%  Usage:
%    1. Edit configuration section below
%    2. Run script
%    3. Results printed to console and displayed in figures
%
%  GitHub: [Add your repository URL here]
% =========================================================================

clear; close all; clc;

%% ========================================================================
%  CONFIGURATION
% =========================================================================
config = create_analysis_config();

% Override specific parameters here if needed:
% config.mouse_tag = '0L2R';
% config.trials_use = 1:18;

%% ========================================================================
%  MAIN ANALYSIS PIPELINE
% =========================================================================

fprintf('\n========================================\n');
fprintf('  ODOR SELECTIVITY ANALYSIS\n');
fprintf('  Mouse: %s\n', config.mouse_tag);
fprintf('  Sessions: %d-%d\n', min(config.sessions_use), max(config.sessions_use));
fprintf('  Selectivity threshold: %.0f%%\n', config.min_fraction_consistent * 100);
fprintf('========================================\n\n');

% Step 1: Analyze neural responses for positive trials
fprintf('Analyzing positive trials...\n');
results_pos = compute_neural_responses(config, 'pos');

% Step 2: Analyze neural responses for negative trials
fprintf('\nAnalyzing negative trials...\n');
results_neg = compute_neural_responses(config, 'neg');

% Step 3: Classify neurons by odor selectivity
fprintf('\nClassifying odor selectivity...\n');
selectivity = classify_odor_selectivity(results_pos, results_neg, config);

% Step 4: Display results
print_selectivity_summary(selectivity, config);

% Step 5: Visualize results
visualize_selectivity_trends(selectivity, config);

fprintf('\n========================================\n');
fprintf('  Analysis Complete!\n');
fprintf('========================================\n\n');


%% ========================================================================
%  CONFIGURATION FUNCTION
% =========================================================================
function config = create_analysis_config()
    % CREATE_ANALYSIS_CONFIG - Creates configuration structure with all parameters
    %
    % Returns:
    %   config - Structure containing all analysis parameters
    
    % Mouse and session parameters
    config.mouse_tag = '0L1R';          % Mouse identifier (e.g., '0L1R', '0L2R', '1L0R')
    config.trials_use = 1:20;           % Trial indices to analyze per session
    config.sessions_use = 1:8;          % Session indices (Day = session + 2)
    
    % Signal processing parameters
    config.fs = 10;                     % Sampling rate (Hz)
    config.baseline_duration = 5;       % Baseline window duration (seconds)
    config.response_duration = 5;       % Response window duration (seconds)
    config.downsample_factor = 10;      % Downsample factor for binning (10 samples → 1)
    
    % Selectivity criteria
    config.min_fraction_consistent = 0.75;  % Minimum fraction of trials for selectivity (0.75 = 75%)
    config.response_threshold_sigma = 2;     % Z-score threshold for response detection
    
    % File naming patterns
    config.behavior_prefix = 'MCA1';    % Behavior file prefix
    config.trace_prefix = 'M1CA1';      % Trace file prefix  
    config.event_prefix = 'M1CA1';      % Event file prefix
    config.data_folder = config.mouse_tag;  % Data folder name
    
    % Derived parameters
    config.n_baseline_samples = config.baseline_duration * config.fs;
    config.n_response_samples = config.response_duration * config.fs;
end


%% ========================================================================
%  NEURAL RESPONSE COMPUTATION
% =========================================================================
function results = compute_neural_responses(config, trial_type)
    % COMPUTE_NEURAL_RESPONSES - Analyzes neural responses for one trial type
    %
    % Inputs:
    %   config      - Configuration structure
    %   trial_type  - 'pos' for positive trials, 'neg' for negative trials
    %
    % Returns:
    %   results - Structure array with per-session results
    
    % Load behavioral data
    behavior = load_behavior_data(config, trial_type);
    
    % Initialize results structure
    results = struct([]);
    
    % Process each session
    for ses = config.sessions_use
        day_number = ses + 2;
        
        % Check if session exists in behavior data
        if ses + 1 > size(behavior.session_blocks, 2)
            warning('Session %d not found in behavior data. Skipping.', ses);
            continue;
        end
        
        % Load neural data
        [traces, events] = load_neural_data(config, day_number);
        if isempty(traces) || isempty(events)
            warning('Session %d: Neural data not found. Skipping.', ses);
            continue;
        end
        
        % Analyze responses for this session
        session_results = analyze_session_responses(traces, events, behavior, ...
                                                    ses, config);
        
        % Store results
        results(ses).day = day_number;
        results(ses).base_mean = session_results.base_mean;
        results(ses).resp_mean = session_results.resp_mean;
        results(ses).resp_gt_base = session_results.resp_gt_base;
        results(ses).valid_mask = session_results.valid_mask;
        
        fprintf('  Session %d (Day %d) analyzed\n', ses, day_number);
    end
end


%% ========================================================================
%  BEHAVIORAL DATA LOADING
% =========================================================================
function behavior = load_behavior_data(config, trial_type)
    % LOAD_BEHAVIOR_DATA - Loads and processes behavioral trial flags
    
    behavior_file = sprintf('%s-%s_behav.csv', config.behavior_prefix, config.mouse_tag);
    
    if ~isfile(behavior_file)
        error('Behavior file not found: %s', behavior_file);
    end
    
    % Read behavioral data
    tab_beh = readtable(behavior_file);
    trial_data = tab_beh{:, :};
    
    % Extract positive and negative trial flags
    trial_pos = trial_data(:, 1);
    trial_neg = trial_data(:, 2);
    trial_pos(isnan(trial_pos)) = 0;
    trial_neg(isnan(trial_neg)) = 0;
    
    % Select appropriate trial type
    switch lower(trial_type)
        case 'pos'
            trial_flags = trial_pos;
        case 'neg'
            trial_flags = trial_neg;
        otherwise
            error('trial_type must be ''pos'' or ''neg''');
    end
    
    % Reshape into sessions
    n_trials_per_session = numel(config.trials_use);
    n_total_trials = numel(trial_flags);
    n_sessions = n_total_trials / n_trials_per_session;
    
    if mod(n_total_trials, n_trials_per_session) ~= 0
        error('Trial count mismatch: expected multiple of %d trials', n_trials_per_session);
    end
    
    % Store in output structure
    behavior.session_blocks = reshape(trial_flags, n_trials_per_session, n_sessions);
    behavior.trial_type = trial_type;
end


%% ========================================================================
%  NEURAL DATA LOADING
% =========================================================================
function [traces, events] = load_neural_data(config, day_number)
    % LOAD_NEURAL_DATA - Loads calcium traces and event times for a session
    
    traces = [];
    events = [];
    
    % Load calcium traces
    trace_file = fullfile(config.data_folder, ...
        sprintf('%s-%s-Day%d-PP-BP-MC-CNMFE-cell-traces.csv', ...
                config.trace_prefix, config.mouse_tag, day_number));
    
    if ~isfile(trace_file)
        warning('Trace file not found: %s', trace_file);
        return;
    end
    
    trace_table = readtable(trace_file);
    traces.time = trace_table{:, 1};
    traces.data = trace_table{:, 2:end};  % Time × Neurons
    traces.data = traces.data - mean(traces.data);  % Center data
    [traces.n_timepoints, traces.n_neurons] = size(traces.data);
    
    % Load event times
    event_file = fullfile(config.data_folder, ...
        sprintf('%s-%s %d.csv', config.event_prefix, config.mouse_tag, day_number));
    
    if ~isfile(event_file)
        warning('Event file not found: %s', event_file);
        traces = [];
        return;
    end
    
    event_table = readtable(event_file);
    
    % Convert event times to sample indices
    n_trials = numel(config.trials_use);
    events.onset_offset = (event_table(1:n_trials, 2:3).Variables - ...
                          event_table{end, 2}) * config.fs;
end


%% ========================================================================
%  SESSION-LEVEL RESPONSE ANALYSIS
% =========================================================================
function results = analyze_session_responses(traces, events, behavior, ses, config)
    % ANALYZE_SESSION_RESPONSES - Computes baseline vs response for each neuron/trial
    
    N = traces.n_neurons;
    K = numel(config.trials_use);
    T = traces.n_timepoints;
    
    % Initialize output arrays
    results.base_mean = nan(N, K);
    results.resp_mean = nan(N, K);
    results.resp_gt_base = false(N, K);
    results.valid_mask = false(N, K);
    
    % Process each neuron and trial
    for neuron = 1:N
        for k = 1:K
            trial = config.trials_use(k);
            
            % Check if this trial is the correct type (pos/neg)
            if behavior.session_blocks(trial, ses + 1) == 0
                continue;
            end
            
            % Get event timing
            on_off = events.onset_offset(trial, :);
            if any(isnan(on_off))
                continue;
            end
            
            onset = round(on_off(1));
            offset = round(on_off(2));
            
            % Define baseline window (before onset)
            base_start = max(1, onset - config.n_baseline_samples + 1);
            base_end = min(T, onset);
            
            % Define response window (before offset)
            resp_start = max(1, offset - config.n_response_samples + 1);
            resp_end = min(T, offset);
            
            % Validate windows
            if base_end <= base_start || resp_end <= resp_start
                continue;
            end
            
            % Extract signal segments
            baseline_signal = traces.data(base_start:base_end, neuron);
            response_signal = traces.data(resp_start:resp_end, neuron);
            
            % Z-score normalization using baseline statistics
            baseline_mean = mean(baseline_signal, 'omitnan');
            baseline_std = std(baseline_signal, 0, 'omitnan');
            
            if baseline_std == 0
                baseline_std = eps;  % Avoid division by zero
            end
            
            z_baseline = (baseline_signal - baseline_mean) / baseline_std;
            z_response = (response_signal - baseline_mean) / baseline_std;
            
            % Downsample (bin) the signals
            n_bins = round(config.baseline_duration * config.fs / config.downsample_factor);
            z_baseline_binned = resample(z_baseline, n_bins, ...
                                        config.baseline_duration * config.fs);
            z_response_binned = resample(z_response, n_bins, ...
                                        config.response_duration * config.fs);
            
            % Compute binned statistics
            base_mean_binned = mean(z_baseline_binned, 'omitnan');
            resp_mean_binned = mean(z_response_binned, 'omitnan');
            base_std_binned = std(z_baseline_binned, 0, 'omitnan');
            
            if ~isfinite(base_mean_binned) || ~isfinite(resp_mean_binned) || ...
               ~isfinite(base_std_binned)
                continue;
            end
            
            if base_std_binned == 0
                base_std_binned = eps;
            end
            
            % Store means
            results.base_mean(neuron, k) = base_mean_binned;
            results.resp_mean(neuron, k) = resp_mean_binned;
            results.valid_mask(neuron, k) = true;
            
            % Detect significant response:
            % Response > baseline if any binned sample exceeds threshold
            for i = 1:length(z_response_binned)
                z_diff = (z_response_binned(i) - z_baseline_binned(i)) / base_std_binned;
                if z_diff > config.response_threshold_sigma
                    results.resp_gt_base(neuron, k) = true;
                    break;
                end
            end
        end
    end
end


%% ========================================================================
%  SELECTIVITY CLASSIFICATION
% =========================================================================
function selectivity = classify_odor_selectivity(results_pos, results_neg, config)
    % CLASSIFY_ODOR_SELECTIVITY - Classifies neurons as pos/neg/both/non-selective
    
    max_sessions = max(config.sessions_use);
    
    % Initialize output structure
    selectivity.pos_only_idx = cell(max_sessions, 1);
    selectivity.neg_only_idx = cell(max_sessions, 1);
    selectivity.both_idx = cell(max_sessions, 1);
    selectivity.none_idx = cell(max_sessions, 1);
    
    selectivity.n_pos = nan(max_sessions, 1);
    selectivity.n_neg = nan(max_sessions, 1);
    selectivity.n_both = nan(max_sessions, 1);
    selectivity.n_none = nan(max_sessions, 1);
    selectivity.n_total = nan(max_sessions, 1);
    
    % Classify neurons for each session
    for ses = config.sessions_use
        % Validate data exists
        if ses > numel(results_pos) || ses > numel(results_neg)
            continue;
        end
        if isempty(results_pos(ses).resp_mean) || isempty(results_neg(ses).resp_mean)
            continue;
        end
        
        % Get response matrices
        pos_responses = results_pos(ses).resp_gt_base;
        pos_valid = results_pos(ses).valid_mask;
        neg_responses = results_neg(ses).resp_gt_base;
        neg_valid = results_neg(ses).valid_mask;
        
        % Check dimensions match
        [N_pos, ~] = size(pos_responses);
        [N_neg, ~] = size(neg_responses);
        
        if N_pos ~= N_neg
            warning('Session %d: Dimension mismatch between pos/neg. Skipping.', ses);
            continue;
        end
        
        N = N_pos;
        selectivity.n_total(ses) = N;
        
        % Calculate response consistency for each neuron
        pos_hits = sum(pos_responses & pos_valid, 2);   % # trials with response
        pos_total = sum(pos_valid, 2);                  % # valid trials
        neg_hits = sum(neg_responses & neg_valid, 2);
        neg_total = sum(neg_valid, 2);
        
        % Compute fractions
        pos_fraction = pos_hits ./ max(pos_total, 1);
        neg_fraction = neg_hits ./ max(neg_total, 1);
        
        % Classify based on consistency threshold
        is_pos_consistent = (pos_total > 0) & (pos_fraction >= config.min_fraction_consistent);
        is_neg_consistent = (neg_total > 0) & (neg_fraction >= config.min_fraction_consistent);
        
        % Find indices for each category
        pos_only = find(is_pos_consistent & ~is_neg_consistent);
        neg_only = find(is_neg_consistent & ~is_pos_consistent);
        both = find(is_pos_consistent & is_neg_consistent);
        none = find(~is_pos_consistent & ~is_neg_consistent);
        
        % Store results
        selectivity.pos_only_idx{ses} = pos_only;
        selectivity.neg_only_idx{ses} = neg_only;
        selectivity.both_idx{ses} = both;
        selectivity.none_idx{ses} = none;
        
        selectivity.n_pos(ses) = numel(pos_only);
        selectivity.n_neg(ses) = numel(neg_only);
        selectivity.n_both(ses) = numel(both);
        selectivity.n_none(ses) = numel(none);
    end
end


%% ========================================================================
%  RESULTS DISPLAY
% =========================================================================
function print_selectivity_summary(selectivity, config)
    % PRINT_SELECTIVITY_SUMMARY - Prints classification results to console
    
    fprintf('\n========================================\n');
    fprintf('  SELECTIVITY CLASSIFICATION RESULTS\n');
    fprintf('========================================\n\n');
    
    for ses = config.sessions_use
        if isnan(selectivity.n_total(ses))
            continue;
        end
        
        N = selectivity.n_total(ses);
        day = ses + 2;
        
        fprintf('Session %d (Day %d) - Total: %d neurons\n', ses, day, N);
        fprintf('  Positive-selective:  %3d cells (%.1f%%)\n', ...
                selectivity.n_pos(ses), 100 * selectivity.n_pos(ses) / N);
        fprintf('  Negative-selective:  %3d cells (%.1f%%)\n', ...
                selectivity.n_neg(ses), 100 * selectivity.n_neg(ses) / N);
        fprintf('  Both-odor response:  %3d cells (%.1f%%)\n', ...
                selectivity.n_both(ses), 100 * selectivity.n_both(ses) / N);
        fprintf('  Non-selective:       %3d cells (%.1f%%)\n\n', ...
                selectivity.n_none(ses), 100 * selectivity.n_none(ses) / N);
    end
end


%% ========================================================================
%  VISUALIZATION
% =========================================================================
function visualize_selectivity_trends(selectivity, config)
    % VISUALIZE_SELECTIVITY_TRENDS - Creates plots of selectivity across sessions
    
    sessions = config.sessions_use;
    
    % Create figure with subplots
    figure('Name', 'Odor Selectivity Trends', 'Color', 'w', 'Position', [100 100 1200 800]);
    
    categories = {'Positive-selective', 'Negative-selective', 'Both-odor', 'Non-selective'};
    data_fields = {'n_pos', 'n_neg', 'n_both', 'n_none'};
    colors = {'#0072BD', '#D95319', '#7E2F8E', '#77AC30'};
    markers = {'o', 's', '^', 'd'};
    
    for i = 1:4
        subplot(2, 2, i);
        data = selectivity.(data_fields{i})(sessions);
        
        plot(sessions, data, [markers{i} '-'], ...
             'LineWidth', 2, 'MarkerSize', 8, ...
             'Color', colors{i}, 'MarkerFaceColor', colors{i});
        
        xlabel('Session Index', 'FontSize', 11, 'FontWeight', 'bold');
        ylabel('Number of Neurons', 'FontSize', 11, 'FontWeight', 'bold');
        title(sprintf('%s Neurons', categories{i}), ...
              'FontSize', 12, 'FontWeight', 'bold');
        grid on;
        box on;
        
        % Set axis properties
        xlim([min(sessions) - 0.5, max(sessions) + 0.5]);
        xticks(sessions);
        set(gca, 'FontSize', 10);
    end
    
    sgtitle(sprintf('Mouse %s - Odor Selectivity Across Sessions', config.mouse_tag), ...
            'FontSize', 14, 'FontWeight', 'bold');
    
    % Create summary figure with all categories
    figure('Name', 'Selectivity Summary', 'Color', 'w', 'Position', [150 150 900 600]);
    hold on;
    
    for i = 1:4
        data = selectivity.(data_fields{i})(sessions);
        plot(sessions, data, [markers{i} '-'], ...
             'LineWidth', 2, 'MarkerSize', 8, ...
             'DisplayName', categories{i}, ...
             'Color', colors{i}, 'MarkerFaceColor', colors{i});
    end
    
    xlabel('Session Index', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('Number of Neurons', 'FontSize', 12, 'FontWeight', 'bold');
    title(sprintf('Mouse %s - All Selectivity Categories', config.mouse_tag), ...
          'FontSize', 14, 'FontWeight', 'bold');
    legend('Location', 'best', 'FontSize', 10);
    grid on;
    box on;
    xlim([min(sessions) - 0.5, max(sessions) + 0.5]);
    xticks(sessions);
    set(gca, 'FontSize', 11);
end

