
function contrast_index = calculate_contrast(I_feature,I_surronding)
% calculates the contrast factor between the brightness of the feature
% (aka inside the river channel) and the brightness outside of the feature
% (aka outside the river channel but within 3.5 the local width)
%   INPUT:
%       I_feature       = brightness of pixels within the channnel
%       I_surronding    = brightness of pixels outside the channel
%   OUTPUT
%       contrast_index  = contrast between inside and outside of channel (from Ajay's draft)

    I_feature = mean(I_feature,'omitnan');
    I_surronding = mean(I_surronding,'omitnan');

    contrast_index = (I_feature - I_surronding)/(I_feature + I_surronding); % after equation 1 in Ajay's draft


end