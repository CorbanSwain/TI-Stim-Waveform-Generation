classdef ModulationSymetry < int8
enumeration
    SIGNAL1 (1)
    SYMETRIC (0)
    SIGNAL2 (-1)
end

methods (Static)
    function isMem = isMember(x)
        try
            ModulationSymetry(x);
        catch 
            isMem = false;
            return
        end
        isMem = true;
    end
end
end