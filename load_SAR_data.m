function [uncorr_sar_file,uncorr_sar_data] = load_SAR_data(sar_name)
% read in data from cub file
%
% INPUT:
%   sar_name = name of pass which will be the first part of the cub file name
%
% OUTPUT:
%   sar_file = filename containing feature data (.cub)
%   sar_data = data from .cub file read in using read_isis external function




uncorr_sar_file = strcat(sar_name,'_reproc1_corramb_bidr.s0nsqt_uncorr.cub');

if ~isfile(fullfile('data','BIDR_SBDR','uncorr',uncorr_sar_file))
    uncorr_sar_file = strcat(sar_name,'_reproc2_corramb_bidr.s0nsqt_uncorr.cub');
end

if ~isfile(fullfile('data','BIDR_SBDR','uncorr',uncorr_sar_file))
    fprintf('...no associated SAR file found for %s\n',sar_name)
    uncorr_sar_file = '';
    uncorr_sar_data = NaN;
    return
end

fprintf('...using %s file for brightness extraction\n',uncorr_sar_file)
uncorr_sar_data = read_isis(uncorr_sar_file); % found on W drive

bad_vals = uncorr_sar_data<-2e38; % no brightness vals ~ -3.4028e+38

if sum(bad_vals(:))
    % fprintf('...SAR Swath %s contains blown out pixels. Removing NaN (~-3e-38) values from brightness pixels.\n',sar_name)
    uncorr_sar_data(bad_vals) = NaN;

end


end