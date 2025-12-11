# Odor Selectivity Analysis

MATLAB pipeline for identifying odor-selective neurons from calcium imaging data during behavioral training sessions.

## Overview

This tool analyzes one-photon calcium imaging data to classify neurons based on their selectivity to positive vs. negative odor cues. It processes multiple training sessions and tracks how selectivity evolves over time.

## Features

- **Automated pipeline**: From raw data to classification with minimal configuration
- **Flexible parameters**: Easy adjustment of selectivity thresholds and time windows
- **Comprehensive output**: Statistical summaries and publication-ready visualizations
- **Modular design**: Well-documented functions for easy customization
- **Session tracking**: Analyze neural selectivity across multiple training days

## Classification Categories

Neurons are classified into four categories based on their response consistency:

1. **Positive-selective**: Respond consistently to positive odors only
2. **Negative-selective**: Respond consistently to negative odors only  
3. **Both-odor responsive**: Respond to both odor types
4. **Non-selective**: No consistent response pattern

## Requirements

- MATLAB R2019b or later
- Signal Processing Toolbox (for `resample` function)

## Data Structure

The pipeline expects the following file organization:

```
project_folder/
├── analyze_odor_selectivity.m     # Main analysis script
├── MCA1-{mouse_tag}_behav.csv     # Behavioral trial flags
└── {mouse_tag}/                    # Data folder for each mouse
    ├── M1CA1-{mouse_tag}-Day3-PP-BP-MC-CNMFE-cell-traces.csv
    ├── M1CA1-{mouse_tag}-Day4-PP-BP-MC-CNMFE-cell-traces.csv
    ├── M1CA1-{mouse_tag} 3.csv     # Event times
    ├── M1CA1-{mouse_tag} 4.csv
    └── ...
```

### Input File Formats

#### Behavioral File (`MCA1-{mouse_tag}_behav.csv`)
- Column 1: Positive trial flags (1 = positive, 0 = not positive)
- Column 2: Negative trial flags (1 = negative, 0 = not negative)
- Rows: All trials from all sessions concatenated

#### Trace Files (`*-cell-traces.csv`)
- Column 1: Time vector (seconds)
- Columns 2+: Fluorescence traces for each neuron
- Rows: Time points

#### Event Files (`* {day}.csv`)
- Columns 1-2: Trial onset and offset times (seconds)
- Last row: Reference time for synchronization

## Quick Start

1. **Clone or download** this repository

2. **Organize your data** according to the structure above

3. **Edit configuration** in `analyze_odor_selectivity.m`:

```matlab
config.mouse_tag = '0L1R';          % Your mouse ID
config.trials_use = 1:20;           % Trials per session
config.sessions_use = 1:8;          % Sessions to analyze
config.min_fraction_consistent = 0.75;  % Selectivity threshold (75%)
```

4. **Run the analysis**:

```matlab
analyze_odor_selectivity
```

## Configuration Parameters

### Basic Parameters
- `mouse_tag`: Mouse identifier (e.g., '0L1R', '0L2R')
- `trials_use`: Trial indices to analyze per session (e.g., 1:20)
- `sessions_use`: Session indices to analyze (e.g., 1:8)

### Signal Processing
- `fs`: Sampling rate in Hz (default: 10)
- `baseline_duration`: Baseline window length in seconds (default: 5)
- `response_duration`: Response window length in seconds (default: 5)
- `downsample_factor`: Binning factor for temporal downsampling (default: 10)

### Selectivity Criteria
- `min_fraction_consistent`: Minimum fraction of trials showing response (default: 0.75 = 75%)
- `response_threshold_sigma`: Z-score threshold for response detection (default: 2)

## Output

### Console Output
- Session-by-session summary with neuron counts and percentages for each category
- Progress indicators during analysis

### Figures
1. **4-panel plot**: Individual trends for each selectivity category
2. **Summary plot**: All categories overlaid for comparison

### Data Structures
- `results_pos`: Per-session responses to positive trials
- `results_neg`: Per-session responses to negative trials
- `selectivity`: Classification results with neuron indices for each category

## Method Details

### Response Detection
1. Extract baseline (before odor onset) and response (before odor offset) windows
2. Z-score normalize using baseline statistics
3. Downsample/bin signals for noise reduction
4. Detect response if `(z_response - z_baseline) / std_baseline > threshold`

### Selectivity Classification
- A neuron is "selective" if it responds in ≥75% of valid trials (configurable)
- Classification based on independent analysis of positive and negative trials

## Customization

### Modify Analysis Pipeline
The code is organized into modular functions:
- `create_analysis_config()`: Set all parameters
- `compute_neural_responses()`: Process calcium traces
- `classify_odor_selectivity()`: Determine selectivity
- `visualize_selectivity_trends()`: Generate plots

### Example: Change Response Criterion
Edit in `analyze_session_responses()`:
```matlab
if z_diff > config.response_threshold_sigma
    results.resp_gt_base(neuron, k) = true;
end
```

### Example: Add New Visualization
Add a new plotting function after `visualize_selectivity_trends()` call in main script.

## Troubleshooting

**"Behavior file not found"**
- Check that `MCA1-{mouse_tag}_behav.csv` exists in the working directory
- Verify `mouse_tag` matches filename exactly

**"Trace file not found"**
- Ensure data folder name matches `mouse_tag`
- Check that day numbers = session + 2 (Session 1 → Day 3)
- Verify file naming format matches expected pattern

**Dimension mismatch warnings**
- Check that same neurons are present in positive and negative trial analyses
- Verify trace files contain consistent neuron counts across sessions

**No neurons classified as selective**
- Try lowering `min_fraction_consistent` threshold
- Check that behavioral flags correctly identify trial types
- Verify calcium signals show visible responses

## Citation

If you use this code in your research, please cite:

```
Nazari, M. (2025). Odor Selectivity Analysis Pipeline. 
GitHub: [Your Repository URL]
```

## Author

**M. Nazari**  
Fall 2025

## License

MIT License - feel free to use and modify for your research.

## Contributing

Contributions welcome! Please open an issue or pull request for:
- Bug fixes
- New features
- Documentation improvements
- Example datasets

## Contact

For questions or issues, please open a GitHub issue or contact [your email].

