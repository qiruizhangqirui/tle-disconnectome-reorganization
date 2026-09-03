function result = calculate_welch_two_sample_statistics(group_one, group_two)
% calculate_welch_two_sample_statistics: Compare two independent groups
%
% Purpose:
% Estimate a two-sided Welch t test, the group one minus group two mean
% difference and confidence interval, and Hedges' g using the pooled sample SD.
%
% Author: Qirui Zhang
% Creation Date: 2026-08-31

% Calculate group descriptives and the group one minus group two difference.
result.N_Group_One = numel(group_one);
result.Mean_Group_One = mean(group_one);
result.SD_Group_One = std(group_one, 0);
result.N_Group_Two = numel(group_two);
result.Mean_Group_Two = mean(group_two);
result.SD_Group_Two = std(group_two, 0);
result.Mean_Difference = result.Mean_Group_One - result.Mean_Group_Two;

% Run the two-sided Welch test and obtain the 95% mean-difference interval.
[~, raw_p, confidence_interval, test_statistics] = ttest2( ...
    group_one, group_two, 'Vartype', 'unequal');

% Standardize the contrast with the pooled sample SD and Hedges correction.
variance_one = var(group_one, 0);
variance_two = var(group_two, 0);
pooled_variance = ( ...
    (result.N_Group_One - 1) * variance_one + ...
    (result.N_Group_Two - 1) * variance_two) / ...
    (result.N_Group_One + result.N_Group_Two - 2);
cohen_d = result.Mean_Difference / sqrt(pooled_variance);
small_sample_correction = 1 - 3 / ( ...
    4 * (result.N_Group_One + result.N_Group_Two) - 9);

result.CI95_Lower = confidence_interval(1);
result.CI95_Upper = confidence_interval(2);
result.Welch_T = test_statistics.tstat;
result.Welch_DF = test_statistics.df;
result.Hedges_g = small_sample_correction * cohen_d;
result.Raw_P = raw_p;
end
