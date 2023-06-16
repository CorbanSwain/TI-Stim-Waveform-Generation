function strOutput = struct2str(s, sName, preSpaces)

VECTOR_LEN_CUTOFF = 32;

switch nargin
   case 1
      sName = 'struct';
      preSpaces = 2;
   case 2
      preSpaces = 2;
end

fnames = fieldnames(s);

if ~isscalar(s)
   strArray = [{sprintf("%s (a %s struct with fields)", sName, ...
      join(string(size(s)), "x"))}; fnames];
   strArray = string(strArray);
else
   nFields = length(fnames);
   maxLen = 0;
   for iField = 1:nFields
      if length(fnames{iField}) > maxLen
         maxLen = length(fnames{iField});
      end
   end
   
   strArray = strings(1, nFields);
   for iField = 1:nFields
      x = s.(fnames{iField});

      if (ischar(x) && isvector(x)) || (isstring(x) && isscalar(x))         
         if ischar(x)
            strArray(iField) = sprintf("'%s'", x);
         else
            strArray(iField) = sprintf('"%s"', x);
         end
      elseif isempty(x)
         if sum(size(x), 'all') == 0
            if isa(x, 'double')
               strArray(iField) = "[]";
            elseif isa(x, 'cell')
               strArray(iField) = "{}";
            else               
               strArray(iField) = sprintf("[] (an empty %s array)", class(x));
            end
         else
            strArray(iField) = sprintf("(a %s empty %s array)", ...
               join(string(size(x)), "-by-"), class(x));
         end
      elseif isstruct(x)
         strArray(iField) = csmu.struct2str(x, '', preSpaces + 4 + maxLen);
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
            strArray(iField) = stringsTemp;
         else
            strArray(iField) = sprintf("[%s]", join(stringsTemp));
            if iscolumn(x)
               strArray(iField) = strcat(strArray(iField), "'");
            end
         end
      else
         if isscalar(x)
            if isenum(x)
               strArray(iField) = sprintf('%s.%s', class(x), string(x));
            else
               strArray(iField) = sprintf("(a %s)", class(x));
            end
         else
            strArray(iField) = sprintf("(a %s %s array)", ...
               join(string(size(x)), "-by-"), class(x));
         end
      end
      strArray(iField) = sprintf("%-*s: %s", maxLen, fnames{iField}, ...
         strArray(iField));
    
   end
   
   strArray = [sName, strArray];
end

strSprintfPrepped = strings(1, length(strArray));
for iStr = 1:length(strArray)
   strSprintfPrepped(iStr) = strrep(strArray(iStr), '\', '\\');
end

strJoined = join(strSprintfPrepped, ['\n', repmat(' ', 1, preSpaces)]);

strOutput = sprintf(strJoined);
end