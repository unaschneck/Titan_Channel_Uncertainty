function [info,fileout]=sar_readsbdr(filename,names_wanted,sflag);

if ~exist('sflag','var');
    sflag = 0;
    outfile = 'None_Saved';
else
    sflag = 1;
    outfile = [filename '.mat'];
end;
%read the appropriate format file
pos = strfind(filename,'sbdr');
pos2 = strfind(filename,'lbdr');
pos3 = strfind(filename,'.TAB');
if length(pos) > 0 | length(pos3) > 0
    type = 'sbdr';
    if exist('SBDR.mat')==2;
        load SBDR.mat
    else;
        fmt = pds_label_parse_v2('SBDR.FMT');
    end;
elseif length(pos2) > 0;
    type = 'lbdr';
    if exist('SBDR.mat')==2;
        load SBDR.mat
    else;
        fmt = pds_label_parse_v2('SBDR.FMT');
    end;
else;
    type = 'unrecognized';
    disp(['Unrecogonized Table Format: ' filename]);
    fileout = nan;
    info = nan;
    return
end;

names = fields(fmt);
% if ~exist('names','var');
%    names = fields(fmt);
% elseif length(names) == 0;
%     names = fields(fmt);
% end;

%read the header
lbl = pds_label_parse_v2(filename);

%get the number of file records and the offset
offset = lbl.recordbytes.*lbl.labelrecords;
switch type
    case{'sbdr'}
        rowbytes = lbl.sbdr_table.rowbytes;
        rows = lbl.sbdr_table.rows;
    case{'lbdr'}
        rowbytes = lbl.lbdr_table.rowbytes;
        rows = lbl.lbdr_table.rows;
    otherwise
        type = 'unrecognized';
        disp(['Unrecogonized Table Format: ' filename]);
        fileout = nan;
        info = nan;
        return
end;

%step through and populate the infor structure
fid = fopen(filename,'r');
if strcmp(names{1},'pds_label_parse_rundate'); 
    temp = names;
    clear names;
    for j=2:length(temp); names{j-1} = temp{j}; end;
end;  
for i=1:length(names);    
    fmt_temp = fmt.(names{i});
    good = 0;
    if exist('names_wanted','var');
        for j=1:length(names_wanted);
            if strcmp(names_wanted{j},fmt_temp.name);
                good = 1;
            end;
        end;
    else;
        good = 1;
    end;
    if good;
        bytes = fmt_temp.bytes;
        switch fmt_temp.datatype;
            case{'pc_real'};
                prec = ['float' num2str(8.*bytes)];
            case{'pc_integer'}
                prec = ['int' num2str(8.*bytes)];
            case{'pc_unsigned_integer'}
                prec = ['uint' num2str(8.*bytes)];                    
            case{'character'};
                prec = 'nan';
            otherwise
                prec = 'nan';;
        end;
        clear val
        if ~strcmp(prec,'nan');
            fseek(fid,offset+fmt_temp.startbyte-1,'bof');
            val = fread(fid,rows,prec,rowbytes-bytes);
            eval(['info.' fmt_temp.name ' = val;']);
        end;
    end;
end;
fclose(fid);
if ~exist('info','var');
    disp(['Could Not Find Any Valid Fields!...']);
    info = nan;
    out_file = nan;
    return;
end;
if sflag; save(outfile,'info'); end;
