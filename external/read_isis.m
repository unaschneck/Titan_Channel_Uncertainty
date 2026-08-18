%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Name: read_isis.m
%
%Description: Program to read Cassini SAR PDS, ISIS, and MER VICAR images
%             (will also read most other PDS images as well)%

%Supp. Formats:    ISIS 2 (BSQ)
%                  ISIS 3 (TILE AND BSQ)
%                  CASSINI SAR BIDR PDS PRODUCTS
%                  MER CAMERA PDS EDRs AND RDRs
%                  MER VICAR PRODUCTS (BSQ)
%
%Inputs:
%    sar_file:     Name of PDS image or ISIS .cub file
%                  NOTE: Both Band Sequential and Tiled Formats are
%                        supported. Hyperspectral cubes work as well
%    samples:      Vector of size [n,1] where the elements are the sample
%                  positions of the regions of interest. May also be
%                  one-dimensional pixel indices of sample array is not
%                  input. If empty, the entire image is read.
%    lines:        Vector of size [n,1] where the elements are the line
%                  positions of the regions of interest. May be left empty
%                  if one-dimensional pixel indices are input. Note that
%                  one-dimensional pixel indices are expected to be in bsq
%                  format. The software will then convert them to
%                  (samples,line,band) in order to locate their position in
%                  the ISIS 3 tiles geometry (if cube format is tile)
%    bands:        Vector of size [n,1] where the elements are the band
%                  positions of the regions of interest. May be left empty.
%    stretch_flag: If set, multiples ISIS cube by core multiplier and adds
%                  core offset. This is useful for stretched cubes.
%                  The default is 1 [set]. Most cubes have 1 for their
%                  multiplier and 0 for their offset.
%    sz_s:         The number of samples in the image
%                  (default is to obtain from header)
%    sz_l:         The number of lines in the image
%                  (default is to obtain from header)
%    sz_b:         The number of bands in the image
%                  (default is to obtain from header)
%                  (sz_b defaults to 1 if no header value is found)
%    offset:       The offset in bytes for the image header
%                  (default is to obtain from header)
%    bytes:        The number of bytes per pixel
%                  (default is to obtain from header)
%    prec:         Data precision [uint8,int16,uint16,float32,etc.]
%                  (default is to obtain from header)
%    byteorder:    Byte order (LSB or MSB)
%                  (default is to obtain from header)
%    format:       Denontes format of data (bsq or tile)
%                  NOTE: bil is currently not supported
%                  (default is to obtain from header)
%
%Outputs:
%        values:  Array of size [numel(samples),1] containing pixel data.
%                 If entire image is read, array is resized to
%                 [sz_s,sz_l,sz_b]
%Function Calls:
%               pds_label_parse_v2.m
%
%Required Data Files:
%                    None
%
%Examples:
%        %%Read in a Band Sequential ISIS3 image%
%         %define a filename
%         filename = 't16_reproc1_corramb_bidr_s0nsqt_corr_crop_isis3.cub';
%         %read the entire file
%         img = sar_readpds_v2(filename);
%
%        %%Obtain pixel value for position (32,48,100) in a hyperspectral
%        %%ISIS 3 cube in Tile format
%         %define a filename
%         filename = 'CM_1537733148_1_ir_isis3.cub';
%         %read pixel (32,48,100)
%         value = sar_readpds_v2(filename,32,48,100);
%
%Change Log:
%           March 2007: Origin Version Created by Alex Hayes
%                       (hayes@gps.caltech.edu)
%           03/04/11:   Header added (Alex Hayes)
%           04/18/11:   Updated to work with ISIS3 tile format (Alex Hayes)
%           04/19/11:   Updated to work with multi-band images (Alex Hayes)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function values = read_isis(sar_file,samples,lines,bands,stretch_flag,...
    sz_s,sz_l,sz_b,offset,bytes,prec,byteorder,format);

%define the stretch flag if needed
%(useful for ISIS cubes which have non-unity core multipliers)
if ~exist('stretch_flag','var') | isempty(stretch_flag);
    stretch_flag = 1;
end;

%get the image sample and line numbers
if ~exist('sz_s','var') | isempty(sz_s) | ...
        ~exist('sz_l','var') | isempty(sz_l)
    
    %load the label data
    info = pds_label_parse_v2(sar_file);
    
    %make sure the label is in a format we can understand
    %CASSINI SAR BIDR PDS FORMAT
    if isfield(info,'image');
        dataformat = 'PDS';
    elseif isfield(info,'qube');
        dataformat = 'ISIS2';
    elseif isfield(info,'isiscube');
        dataformat = 'ISIS3';
    elseif isfield(info,'fileformat') & strcmp(info.fileformat,'vicar');
        dataformat = 'VICAR';
    else
        disp('Data Format Not Supported (Not PDS, ISIS2, or ISIS3)!');
        disp(['Please run again with manual input all data format '...
            'paramters (see header)']);
        disp('Returning...');
        values = -1;
        return;
    end;
    
    %get the sample, line, and band numbers from the label
    %CASSINI SAR BIDR PDS FORMAT
    if isfield(info,'image');
        sz_s = info.image.linesamples;
        sz_l = info.image.lines;
        if isfield(info.image,'bands');
            sz_b = info.image.bands;
        else
            sz_b = 1; %Number of bands not nominally listed in the header
        end;
        %VICAR
    elseif isfield(info,'fileformat') & strcmp(info.fileformat,'vicar');
        sz_s = info.ns;
        sz_l = info.nl;
        sz_b = info.nb;
        %ISIS 3
    elseif isfield(info,'isiscube');
        sz_s = info.isiscube.core.samples;
        sz_l = info.isiscube.core.lines;
        sz_b = info.isiscube.core.bands;
        if iscell(sz_b);
            temp = sz_b;
            clear sz_b;
            sz_b = str2num(temp{1});
        end;
        %ISIS 2
    else;
        sz = info.qube.coreitems;
        if ischar(sz);
            sz = str2num(info.qube.coreitems);
        end;
        sz_s = sz(1);
        sz_l = sz(2);
        sz_b = sz(3);
    end;
end

%make sure sz_b was input (set to 1 if not)
if ~exist('sz_b','var') | isempty(sz_b);
    disp('sz_b not input, setting default of sz_b=1');
    sz_b = 1;
end;

%calculate the header offset
if ~exist('offset','var') | isempty(offset);
    if ~exist('info','var') | isempty(info);
        info = pds_label_parse_v2(sar_file);
    end;
    %ISIS 2
    if isfield(info,'cqube');
        offset = info.recordbytes.*(info.cqube-1);
        %ISIS 3
    elseif isfield(info,'isiscube');
        offset = info.isiscube.core.startbyte-1;
        %VICAR
    elseif isfield(info,'fileformat') & strcmp(info.fileformat,'vicar');
        offset = info.lblsize;
        %CASSINI SAR BIDR PDS FORMAT
    elseif isfield(info,'image');
        offset = info.recordbytes.*info.labelrecords;
    else
        disp('Data Format Not Supported (Not PDS, ISIS2, or ISIS3)!');
        disp(['Please run again with manual input all data format '...
            'paramters (see header)']);
        disp('Returning...');
        values = -1;
        return;
    end;
end;

%get the number of bytes
if ~exist('bytes','var') | isempty(bytes);
    if ~exist('info','var') | isempty(info);
        info = pds_label_parse_v2(sar_file);
    end;
    %get bytes and precision
    %VICAR
    if isfield(info,'fileformat') & strcmp(info.fileformat,'vicar');
        pixformat = info.format;
        switch pixformat
            case{'real'}
                bytes =4;
        end;        
    %ISIS 2
    elseif isfield(info,'qube');
        bytes = info.qube.coreitembytes;
        %ISIS 3
    elseif isfield(info,'isiscube');
        type = info.isiscube.core.type;
        switch type
            case{'unsignedbyte'}
                bytes = 1;
            case{'signedword'}
                bytes = 2;
            case{'real'}
                bytes = 4;
        end;
        %CASSINI SAR BIDR PDS FORMAT
    elseif isfield(info,'image');
        bytes = info.image.samplebits./8;
    else
        disp('Data Format Not Supported (Not PDS, ISIS2, or ISIS3)!');
        disp(['Please run again with manual input all data format '...
            'paramters (see header)']);
        disp('Returning...');
        values = -1;
        return;
    end;    
end;

%get the data precision
if ~exist('prec','var') | isempty(prec);
    if ~exist('info','var') | isempty(info);
        info = pds_label_parse_v2(sar_file);
    end;    
    %ISIS 2
    if isfield(info,'qube');
        sampletype = info.qube.coreitemtype;
        %ISIS 3
    elseif isfield(info,'isiscube');
        sampletype = info.isiscube.core.type;
        %VICAR
    elseif isfield(info,'fileformat') & strcmp(info.fileformat,'vicar');
        sampletype = info.format;
        %CASSINI SAR PDS FORMAT
    elseif isfield(info,'image');
        sampletype = info.image.sampletype;
    else
        disp('Data Format Not Supported (Not PDS, ISIS2, or ISIS3)!');
        disp(['Please run again with manual input all data format '...
            'paramters (see header)']);
        disp('Returning...');
        values = -1;
        return;
    end;
    switch sampletype;
        case{'real'}
            prec = 'float32';
            byteorder = 'lsb';
        case{'ieee_real'}
            prec = 'float32';
            byteorder = 'msb';
        case{'pc_real'}
            prec = 'float32';
        case{'mac_real'}
            prec = 'float32';
            byteorder ='msb';
        case{'sun_real'}
            byteorder = 'msb';
            prec = 'float32';
        case{'msb_integer'};
            prec = ['int' num2str(info.image.samplebits)];
            byteorder = 'msb';
        case{'signedword'};
            prec = 'int16';
        case{'unsignedbyte'}
            prec = 'uint8';
        case{'integer'}
            prec = ['int' num2str(info.image.samplebits)];
        case{'unsigned integer'}
            prec = ['uint' num2str(info.image.samplebits)];
        case{'pc_integer'}
            prec = ['int16'];
            stretch_flag = 1;
        otherwise
            disp(['Invalid Sample Type: ' sampletype]);
            disp('Returning...');
            values = -1;
    end;
end;

%get the byte order
if exist('byteorder','var') & ~isempty(byteorder);
    byteorder = lower(byteorder);
    switch byteorder;
        case{'lsb'}
            endian = 'ieee-le';
        otherwise;
            endian = 'ieee-be';
    end;
else;
    if ~exist('info','var') | isempty(info);
        info = pds_label_parse_v2(sar_file);
    end;
    %ISIS3
    if isfield(info,'isiscube');
        byteorder = info.isiscube.core.byteorder;
        %ISIS2 and PDS DATA DEFAUL TO LSB??? (NOT IN HEADER!!)
        %ISIS2
    elseif isfield(info,'qube');
        byteorder = 'lsb';
        %CASSINI SAR PDS
    elseif isfield(info,'image');
        byteorder = 'lsb';
    else
        disp('Data Format Not Supported (Not PDS, ISIS2, or ISIS3)!');
        disp(['Please run again with manual input all data format '...
            'paramters (see header)']);
        disp('Returning...');
        values = -1;
        return;
    end;
    byteorder = lower(byteorder);
    switch byteorder;
        case{'lsb'}
            endian = 'ieee-le';
        otherwise;
            endian = 'ieee-be';
    end;
end;

%open the file
fid = fopen(sar_file,'r',endian);

%determine if we are reading the entire file
if (~exist('samples','var') | isempty(samples)) & ...
        (~exist('lines','var')   | isempty(lines)  );
    entire_file = 1;
else;
    entire_file = 0;
end;

%determine if we are in tile or bsq format
%default is to use bsq in absence of info structure)
if ~exist('format','var') | isempty(format);
    format = 'bsq';
end;
if ~exist('info','var') | isempty(info);
    format = 'bsq';
else;
    %ISIS 3
    if isfield(info,'isiscube');
        format = info.isiscube.core.format;
        if strcmp(format,'bandsequential');
            format = 'bsq';
        end;
    end;
end;

%if we are using the tile format, get the tile dimensions
if strcmp(format,'tile');
    tile_s = info.isiscube.core.tilesamples;
    tile_l = info.isiscube.core.tilelines;
end;

%read the data;
if entire_file
    fseek(fid,offset,'bof');
    
    %reshape the array into the expected size
    switch format
        case {'bsq'};
            values = fread(fid,sz_s.*sz_l*sz_b,prec);
            values = reshape(values,sz_s,sz_l,sz_b);
            disp([endian ' ' prec ' ' num2str(offset)]);
        case {'tile'};
            %we need to order the values based on the tile size
            
            %define the number of tiles in samples or lines
            n_tiles_s = ceil( sz_s ./ tile_s );
            n_tiles_l = ceil( sz_l ./ tile_l );
            n_tiles = n_tiles_s .* n_tiles_l;
            %define the size of the edge tiles (not needed, files zero padded!)
            %tile_s_edge = sz_s - tile_s.*(n_tiles_s -1);
            %tile_l_edge = sz_l - tile_l.*(n_tiles_l - 1);
            %disp([tile_s_edge,tile_l_edge]);
            
            %create the zero padded output array
            values = zeros(n_tiles_s.*tile_s,n_tiles_l.*tile_l,sz_b);
            
            %each (tile_s*tile_l) elements are stored in a band sequential
            %format
            for k=1:sz_b;
                for j=1:n_tiles_l;
                    for i=1:n_tiles_s;
                        %%get the number of columns in the tile
                        %if i==n_tiles_s;
                        %%if i==1;
                        %     len_s = tile_s;%tile_s_edge;
                        %else;
                        %     len_s = tile_s;
                        %end;
                        %get the number of rows in the tile
                        if j==n_tiles_l;
                            %if j==1;
                            len_l = tile_l;%tile_l_edge;
                        else;
                            len_l = tile_l;
                        end;
                        temp = fread(fid,tile_s.*tile_l,prec);
                        temp = reshape(temp,tile_s,tile_l);
                        %place the tile into the array
                        values((1+(tile_s.*(i - 1))):( ((tile_s.*(i-1)+tile_s))),...
                            1+(tile_l.*(j - 1)):((tile_l.*(j-1)+tile_l)),k ) = ...
                            temp;
                    end;
                end;
            end;
            %the image is zero-padded to keep tile size constant, remove zeros
            values_old = values;
            values = zeros(sz_s,sz_l,sz_b);
            for k=1:sz_b;
                values(:,:,k) = values_old(1:sz_s,1:sz_l,k);
            end;
            clear values_old;
    end;
else;
    switch format
        case {'bsq'};
            %determine if we are inputting samples/line or straight index
            if exist('samples','var') & ~exist('lines','var');
                inds = samples - 1;
                clear samples;
            else;
                %round the samples and lines to the nearest integer
                samples = round(samples);
                lines = round(lines);
                %if the bands variable exists, make sure its an integer too
                if exist('bands','var') & ~isempty(bands);
                    bands = round(bands);
                    %convert to indices
                    inds = sub2ind([sz_s,sz_l,sz_b],samples,lines,bands)-1;
                    clear samples lines bands
                else;
                    %convert to indices
                    inds = sub2ind([sz_s,sz_l],samples,lines)-1;
                    clear samples lines
                end;
            end;
            %sort the indices
            [inds,sind] = sort(inds);
            %create the reverse transform
            [trash,sind2] = sort(sind);
            %find the contiguous blocks
            diffvals =diff(inds);
            diffvals(2:(end+1)) = diffvals;
            diffvals(1) = 0;
            inds_diff = find(diffvals ~= 1);
            szout = length(inds);
            inds_diff(end+1) = szout+1;
            %create the output array
            values = zeros(1,szout);
            %step through and read the data in contiguous blocks
            for i=1:(length(inds_diff)-1);
                %move to position
                fseek(fid,offset+inds(inds_diff(i)).*bytes,'bof');
                %read the data
                len = inds_diff(i+1)-inds_diff(i);
                temp = fread(fid,len,prec);
                %fill in array
                values(inds_diff(i):(inds_diff(i+1)-1)) = temp;
            end;
        case{'tile'}
            %if one component indices are input, convert to
            %sample,line,band
            if exist('samples','var') & ~exist('lines','var');
                inds = samples - 1;
                [samples,lines,bands] = ind2sub([sz_s,sz_l,sz_b],inds+1);
                clear samples;
            end;
            %make sure both sample and line are defined
            if ~exist('samples','var') | isempty(samples) | ...
                    ~exist('lines','var')   | isempty(lines);
                disp('Must Define Sample and Line in order to use Tile Format!');
                disp('Returning...');
                values = -1;
                return;
            else;
                %set the bands array to 1 if it does not exist
                if ~exist('bands','var') | isempty(bands);
                    bands = 1;
                end;
                
                %determine the total tile size
                tile_s_tot = ceil(sz_s./tile_s);
                tile_l_tot = ceil(sz_l./tile_l);
                tile_s_num = ceil(samples./tile_s);
                tile_l_num = ceil(lines./tile_l);
                
                %determine the starting positions for each of the starting tiles
                s_start = 1+(tile_s.*(tile_s_num - 1));
                l_start = 1+(tile_l.*(tile_l_num - 1));
                ind_start = tile_s .*  tile_l .* tile_s_tot .* (tile_l_num-1) + ...
                    tile_s.*tile_l .* (tile_s_num - 1) + ...
                    tile_s_tot.*tile_l_tot.*tile_s.*tile_l.*(bands-1);
                
                %determine the positions within each tile
                ind_tile = sub2ind([tile_s,tile_l],samples-s_start+1,lines-l_start+1);
                %determine the indices within the file
                inds = ind_start + ind_tile -1;
                clear samples lines
                %sort the indices
                [inds,sind] = sort(inds);
                %create the reverse transform
                [trash,sind2] = sort(sind);
                %find the contiguous blocks
                diffvals =diff(inds);
                diffvals(2:(end+1)) = diffvals;
                diffvals(1) = 0;
                inds_diff = find(diffvals ~= 1);
                szout = length(inds);
                inds_diff(end+1) = szout+1;
                %create the output array
                values = zeros(1,szout);
                %step through and read the data in contiguous blocks
                for i=1:(length(inds_diff)-1);
                    %move to position
                    fseek(fid,offset+inds(inds_diff(i)).*bytes,'bof');
                    %read the data
                    len = inds_diff(i+1)-inds_diff(i);
                    temp = fread(fid,len,prec);
                    %fill in array
                    values(inds_diff(i):(inds_diff(i+1)-1)) = temp;
                end;
            end;
    end;
    
end;
%resort the values
if exist('sind2','var');
    values = values(sind2);
end;
%stretch the data values if neccessary
if stretch_flag;
    if ~exist('info','var') | isempty(info);
        info = pds_label_parse_v2(sar_file);
    end;
    %ISIS3
    if isfield(info,'isiscube');
        values = values.*info.isiscube.core.multiplier + ...
            info.isiscube.core.base;
        %ISIS 2
    elseif isfield(info,'qube');
        values = values.*info.qube.coremultiplier + info.qube.corebase;
        %CASSINI SAR BIDR PDS
    elseif isfield(info,'image');
        values = values;
        %VICAR
    elseif isfield(info,'fileformat') & strcmp(info.fileformat,'vicar');
        values = values;
    else
        disp('Data Format Not Supported (Not PDS, ISIS2, or ISIS3)!');
        disp(['Was Not Able to Stretch Result using multiplier and offset']);
        disp('Returning...');        
    end;
end;
%close the file
fclose(fid);
end