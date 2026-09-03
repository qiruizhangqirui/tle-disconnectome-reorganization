function [lower_whisker, upper_whisker] = ...
        calculate_box_whisker_limits(values)
%CALCULATE_BOX_WHISKER_LIMITS Return Tukey boxplot whisker endpoints.

values = values(isfinite(values));
if isempty(values)
    error('Block1:EmptyNetworkValues', ...
        'A network has no finite values.');
end

quartiles = quantile(values, [0.25, 0.75]);
interquartile_range = quartiles(2) - quartiles(1);
lower_fence = quartiles(1) - 1.5 * interquartile_range;
upper_fence = quartiles(2) + 1.5 * interquartile_range;
nonoutlier_values = values( ...
    values >= lower_fence & values <= upper_fence);

lower_whisker = min(nonoutlier_values);
upper_whisker = max(nonoutlier_values);
end
