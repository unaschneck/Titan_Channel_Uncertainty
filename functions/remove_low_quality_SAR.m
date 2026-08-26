function [filtered_names,filtered_segs] = remove_low_quality_SAR(pos_swaths,pos_seg)
% Only look at data points for which there is High Res SAR data (radar modes 3 and 11)
%
% INPUT:
%   pos_swaths     = swath numbers that correspond to an observation (e.g., 55, 56)
%   pos_seg        = associated segment number for swath ID
%
% OUTPUT:
%   filtered_names = swath numbers that correspond to an observation AND have high res
%                    radar data available (e.g., 56)
%   filtered_segs  = associated segment numbers that correspond to filtered_names


     
     filtered_names = [];
     filtered_segs  = [];
        for i =  1:numel(pos_swaths)
        
            sar_name = "t" + string(pos_swaths(i));
            
            sar_file_name_v1 = strcat(sar_name,'_reproc1_allstuffed.sbdr');
            sar_file_name_v2 = strcat(sar_name,'_reproc2_allstuffed.sbdr');

            if ~isfile(fullfile('data','BIDR_SBDR','allstuffed',sar_file_name_v2))
            
                if ~isfile(fullfile('data','BIDR_SBDR','allstuffed',sar_file_name_v1))
                    fprintf('No data available for %s\n',sar_name);
                    continue              
                else
                    sar_file_name = sar_file_name_v1;
                    
                end
            else
                sar_file_name = sar_file_name_v2;
            end
    

            fprintf('...file found for %s: %s\n',sar_name,sar_file_name)
            info = sar_readsbdr(char(sar_file_name));
            radar_modes = unique(info.radar_mode);
            if ismember(3,radar_modes) || ismember(11,radar_modes)
                filtered_names(end+1) = pos_swaths(i);
                filtered_segs(end+1) =  pos_seg(i);
            else
                fprintf('...pass %s segment %i removed from list because no high quality SAR data availabilite\n',sar_name,pos_seg(i));
            end
        end

        % cd(cwd)

end
        
        




