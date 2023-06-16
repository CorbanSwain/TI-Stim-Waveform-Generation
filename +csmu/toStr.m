function s = toStr(x)
VECTOR_LEN_CUTOFF = 10;

if (ischar(x) && isvector(x)) || (isstring(x) && isscalar(x))
    if ischar(x)
        s = sprintf("'%s'", x);
    else
        s = sprintf('"%s"', x);
    end
elseif isempty(x)
    if sum(size(x), 'all') == 0
        if isa(x, 'double')
            s = "[]";
        elseif isa(x, 'cell')
            s = "{}";
        else
            s = sprintf("[] (an empty %s array)", class(x));
        end
    else
        s = sprintf("(a %s empty %s array)", ...
            join(string(size(x)), "-by-"), class(x));
    end
elseif (isnumeric(x) || islogical(x)) && isvector(x)
    if length(x) > (VECTOR_LEN_CUTOFF + 1)
        halfCut = floor(VECTOR_LEN_CUTOFF / 2);
        xTemp = x([1:halfCut, (end - halfCut + 1):end]);
        stringsTemp = string(xTemp(:)');
        stringsTemp = [stringsTemp(1:halfCut), "...", ...
            stringsTemp((halfCut+1):end)];
    else
        stringsTemp = string(x(:)');
    end
    stringsTemp(ismissing(stringsTemp)) = "NaN";
    if isscalar(stringsTemp)
        s = stringsTemp;
    else
        s = sprintf("[%s]", join(stringsTemp));
        if iscolumn(x)
            s = strcat(s, "'");
        end
    end
elseif iscell(x) && isvector(x) && (length(x) < VECTOR_LEN_CUTOFF)
    s = '{';
    if iscolumn(x)
        joinStr = ';';
    else
        joinStr = ',';
    end
    for iEle = 1:length(x)
        if iEle == 1
            commaStr = '';
        else
            commaStr = [joinStr, ' '];
        end
        s = [s, commaStr, char(csmu.toStr(x{iEle}))];
    end
    s = [s, '}'];
else
    if isscalar(x)
        if isenum(x)
            s = sprintf('%s.%s', class(x), string(x));
        else
            s = sprintf("(a %s)", class(x));
        end
    else
        s = sprintf("(a %s %s array)", ...
            join(string(size(x)), "-by-"), class(x));
    end
end
s = char(s);
end