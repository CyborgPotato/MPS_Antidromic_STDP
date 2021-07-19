classdef HH_Modified
    %HH_MODIFIED Modified Multi-compartment Hodgkin Huxley Model
    %   Multicompartment HH model for simulating antidromic STDP in
    %   motorneuron pools
    
    properties
        V % Membrane voltage
        
    end
    
    methods
        function obj = HH_Modified(inputArg1,inputArg2)
            %HH_MODIFIED Construct an instance of this class
            %   Detailed explanation goes here
            obj.Property1 = inputArg1 + inputArg2;
        end
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end

