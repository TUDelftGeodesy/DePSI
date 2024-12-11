function data = freadcpx(fid,Nrows,Ncols,data_format);

% function to read a complex file (in buffers)
% 
% data_format can be 'int16','float32','cpxint16' or 'cpxfloat32'
%
% In the future, a number of checks on the input should be added

if strcmp('int16',data_format)|strcmp('float32',data_format)
elseif strcmp('cpxint16',data_format)
  data_format = 'int16';
elseif strcmp('cpxfloat32',data_format)
  data_format = 'float32';
else
  error('A wrong data format was specified');
end


[data,count] = fread(fid,[2*Ncols Nrows],data_format);
data = complex(data(1:2:2*Ncols,:),data(2:2:2*Ncols,:)).';
% dot is important because of non-conjugate transpose
