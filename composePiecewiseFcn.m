function f = composePiecewiseFcn(functions, values)
    function output = pwf(x, functions, values)        
        output = zeros(size(x), 'like', x);
        output(:) = NaN;
        matched = false(size(x));
        for iVal = 1:length(values)
            query = (x <= values(iVal)) & ~matched;               
            if any(query, 'all')
                matched = matched | query;
                output(query) = functions{iVal}(x(query));
            end
            if all(matched, 'all')
                break
            end
        end
    end
f = @(x) pwf(x, functions, values);
end