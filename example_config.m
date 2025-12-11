%% ========================================================================
%  EXAMPLE CONFIGURATION FILE
%  
%  Copy the relevant section below into analyze_odor_selectivity.m
%  to customize parameters for your analysis.
% =========================================================================

%% EXAMPLE 1: Standard Analysis (0L1R mouse)
config.mouse_tag = '0L1R';
config.trials_use = 1:20;
config.sessions_use = 1:8;
config.min_fraction_consistent = 0.75;  % 75% threshold
config.response_threshold_sigma = 2;     % 2-sigma criterion

%% EXAMPLE 2: Mouse with fewer trials (0L2R)
config.mouse_tag = '0L2R';
config.trials_use = 1:18;                % Only 18 trials per session
config.sessions_use = 1:8;
config.min_fraction_consistent = 0.75;
config.response_threshold_sigma = 2;

%% EXAMPLE 3: More stringent selectivity (100% of trials)
config.mouse_tag = '1L0R';
config.trials_use = 1:20;
config.sessions_use = 1:8;
config.min_fraction_consistent = 1.0;    % Require ALL trials
config.response_threshold_sigma = 2;

%% EXAMPLE 4: More lenient selectivity (50% of trials)
config.mouse_tag = '0L1R';
config.trials_use = 1:20;
config.sessions_use = 1:8;
config.min_fraction_consistent = 0.5;    % Only 50% of trials
config.response_threshold_sigma = 2;

%% EXAMPLE 5: Higher response threshold (3-sigma)
config.mouse_tag = '0L1R';
config.trials_use = 1:20;
config.sessions_use = 1:8;
config.min_fraction_consistent = 0.75;
config.response_threshold_sigma = 3;     % More conservative detection

%% EXAMPLE 6: Analyze subset of sessions (early training)
config.mouse_tag = '0L1R';
config.trials_use = 1:20;
config.sessions_use = 1:4;               % Only first 4 sessions
config.min_fraction_consistent = 0.75;
config.response_threshold_sigma = 2;

%% EXAMPLE 7: Analyze subset of sessions (late training)
config.mouse_tag = '0L1R';
config.trials_use = 1:20;
config.sessions_use = 5:8;               % Only last 4 sessions
config.min_fraction_consistent = 0.75;
config.response_threshold_sigma = 2;

%% EXAMPLE 8: Custom signal processing parameters
config.mouse_tag = '0L1R';
config.trials_use = 1:20;
config.sessions_use = 1:8;
config.fs = 10;                          % Sampling rate
config.baseline_duration = 3;            % Shorter baseline (3s instead of 5s)
config.response_duration = 3;            % Shorter response window
config.downsample_factor = 5;            % Less aggressive binning
config.min_fraction_consistent = 0.75;
config.response_threshold_sigma = 2;

%% ========================================================================
%  HOW TO USE THESE EXAMPLES
% =========================================================================
%
% 1. Choose the example that best matches your analysis needs
% 2. Copy the relevant config lines
% 3. Paste into analyze_odor_selectivity.m after the line:
%       config = create_analysis_config();
% 4. The copied lines will override the default values
%
% Example usage in analyze_odor_selectivity.m:
%
%   config = create_analysis_config();
%   
%   % Override for mouse 0L2R with 18 trials
%   config.mouse_tag = '0L2R';
%   config.trials_use = 1:18;
%
% =========================================================================

