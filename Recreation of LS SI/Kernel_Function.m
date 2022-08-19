classdef Kernel_Function
    properties
        Freq {mustBeNumeric}
        Theta 
    end
    methods
        function obj = Kernel_Function(omega, theta1, theta2)
            obj.Freq = omega;
            obj.Theta = [theta1; theta2];
        end
        function k = calculateKernel(obj, arg_1, arg_2)
            k = zeros(length(arg_1),length(arg_2));
            for i = 1:length(arg_1)
                for j = 1:length(arg_2)
                    k(i,j) = exp(-1*(sin(pi()*(arg_1(i)-arg_2(j)))/obj.Freq)^2/(2*obj.Theta(2)))*exp(-1*(arg_1(i)-arg_2(j))^2/(2*obj.Theta(1)));
                end
            end
        end
    end
end