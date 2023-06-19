function f = composePiecewiseFcn(functions, values)
    function output = pwf(x, functions, values)        
        output = zeros(size(x), 'like', x);
        output(:) = NaN;
        
        if isempty(x)
            return
        end

        matched = false(size(x));
        for iVal = 1:length(values)
            query = (x <= values(iVal)) & ~matched;               
            if ~any(query, 'all')
                continue
            end
            matched = matched | query;
            output(query) = functions{iVal}(x(query));
            if all(matched, 'all')
                return
            end
        end
    end

f = @(x) pwf(x, functions, values);
end