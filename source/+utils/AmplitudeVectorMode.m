classdef AmplitudeVectorMode < int8
enumeration
    EXACT (1)
    REPEAT (2)
end

methods (Static)
    function isMem = isMember(x)
        try
            utils.AmplitudeVectorMode(x);
        catch 
            isMem = false;
            return
        end
        isMem = true;
    end
end
end