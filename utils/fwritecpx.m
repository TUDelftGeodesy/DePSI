function fwritecpx(fid,data,data_format)

% function to write a complex file (in buffers)
% 
% data_format can be 'int16','float32','cpxint16' or 'cpxfloat32'
%
% In the future, a number of checks on the input should be added
%

if strcmp('int16',data_format)|strcmp('float32',data_format)
elseif strcmp('cpxint16',data_format)
  data_format = 'int16';
elseif strcmp('cpxfloat32',data_format)
  data_format = 'float32';
else
  error('A wrong data format was specified');
end

data = data.';
data=[real(data), imag(data)];
data=reshape(data,prod(size(data))/2,2).';

fwrite(fid,data,data_format);
