function [pre_bold, post_bold, pre_gmv, post_gmv] = ...
    select_regional_features(feature_data, orientation_index)
% Purpose:
% Select the four regional feature arrays for one approved atlas
% orientation.
%
% Author: Qirui Zhang
% Created: 1 September 2026
%
% feature_data contains the approved PRE and POST BOLD and GMV arrays for
% Original and Flipped_LR atlas orientations.
% orientation_index uses 1 for Original and 2 for Flipped_LR.
%
% pre_bold and post_bold are time points by 116 ROI arrays.
% pre_gmv and post_gmv are 116 by 1 regional GMV arrays in mm3.
% This function selects named arrays without transformation or missing value
% omission.

% Select the four arrays for the approved atlas orientation.
if orientation_index == 1
    pre_bold = feature_data.bold_pre_original;
    post_bold = feature_data.bold_post_original;
    pre_gmv = feature_data.gmv_pre_original;
    post_gmv = feature_data.gmv_post_original;
else
    pre_bold = feature_data.bold_pre_flipped_lr;
    post_bold = feature_data.bold_post_flipped_lr;
    pre_gmv = feature_data.gmv_pre_flipped_lr;
    post_gmv = feature_data.gmv_post_flipped_lr;
end
end
