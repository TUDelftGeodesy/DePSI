function value = getparameter(param,resdata,skip,nonum_flag)
    
    if nargin < 4;
        nonum_flag='';
    end
    
    ind = strmatch(param,resdata);
    
    if isempty(ind);
        error('EXIT: parameter [%s] not in DORIS result file',param);
    elseif isempty(nonum_flag);
        value = str2double(char(resdata(ind(end)+skip)));
    else
        value = char(resdata(ind(end)+skip));
    end
    
end

