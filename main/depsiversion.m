function ver = depsiversion()
% returns the last change SVN revision and date of depsi folder
curfile = mfilename('fullpath');
sep = strfind(curfile, filesep);
mainpath = curfile(1:sep(end));
[status,info] = system(['svn info ' mainpath]);
assert(status == 0, ['svn command "svn info ' mainpath '" failed']);
info = fields_colon(info);
ver = ['DePSI svn rev. ' info.LastChangedRev ', ' info.LastChangedDate];
end

function fields = fields_colon( in )
fields = [];
% tokens = regexp( sprintf('\n%s\n', in), ... % split block.data in two fields:
%   [ '\n(?<key>[^\n:]+):' ... % #1=key, excluding the terminating colon
%     '[ \t]*(?<value>[^\n]*[^:]*)' ... % #2=value, leading whitespace stripped,
%     '(?=\n[^\n:]+:|\n[*])' ], ... % plus adjacent non key:value lines.
%     'names' );
tokens = regexp( sprintf('\n%s\n', in), ... % split block.data in two fields:
  [ '\n(?<key>[^\n:]+):' ... % #1=key, excluding the terminating colon
    '[ \t]*(?<value>[^\n]*)' ], ... % #2=value, leading whitespace stripped
    'names' );
for token = tokens
  fields = setfield( fields, regexprep( token.key, '[^a-zA-Z0-9]+', '' ), token.value );
end
end