function [uncorr_sar_file,uncorr_sar_data] = load_SAR_data(sar_name,seg_number)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read in data from cub file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT:
%   sar_name          = name of pass which will be the first part of the cub file name
%   seg_number        = segment number for pass 
%
% OUTPUT:
%   uncorr_sar_file   = filename containing feature data (.cub)
%   uncorr_sar_data   = data from .cub file read in using read_isis external function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if seg_number == 1 % the main segment is not always the data for the river, check the segment number if its not listed as 1
    uncorr_sar_file = strcat(sar_name,'_reproc1_corramb_bidr.s0nsqt_uncorr.cub');
    
    if ~isfile(fullfile('data','BIDR_SBDR','uncorr',uncorr_sar_file))
        uncorr_sar_file = strcat(sar_name,'_reproc2_corramb_bidr.s0nsqt_uncorr.cub');
    end

else
    seg_file = dir(fullfile('data','BIDR_SBDR','uncorr', [sar_name, '*','seg',num2str(seg_number),'*'])); % finds associated file that includes tXX and segX in the name
    if ~isempty(seg_file)
        uncorr_sar_file = seg_file.name;
    else
        fprintf('...no associated sar file found for %s segment %i\n',sar_name,seg_number)
        uncorr_sar_file = '';
        uncorr_sar_data = NaN;
        return

    end

end

if ~isfile(fullfile('data','BIDR_SBDR','uncorr',uncorr_sar_file))
        fprintf('...no associated sar file found for %s segment %i\n',sar_name,seg_number)
        uncorr_sar_file = '';
        uncorr_sar_data = NaN;
        return
end

fprintf('...using %s file for brightness extraction\n',uncorr_sar_file)
uncorr_sar_data = read_isis(uncorr_sar_file);                              % found on W drive

bad_vals = uncorr_sar_data<-2e38; % no brightness vals ~ -3.4028e+38

if sum(bad_vals(:))
    % fprintf('...SAR Swath %s contains blown out pixels. Removing NaN (~-3e-38) values from brightness pixels.\n',sar_name)
    uncorr_sar_data(bad_vals) = NaN;

end


end