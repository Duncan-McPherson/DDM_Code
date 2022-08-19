classdef GPR
    properties
        x_data
        y_data
        covar
        inv_covar_input
        noise 
        memory
    end
    methods
        function obj = GPR(X_Train_Data, Y_Train_Data, kernel, Noise)
            %All inputs are in 1 by n matrixs
            obj.x_data = X_Train_Data;
            obj.y_data = Y_Train_Data;
            obj.covar = kernel;
            obj.noise = Noise;

            %Calculate other matrixs
            obj.inv_covar_input = inv(kernel.calculateKernel(X_Train_Data, X_Train_Data) + (Noise+3*10^-7)*eye(size(X_Train_Data,1)));
        end

        function [avg, obj] = Prediction(obj, test_Values)
            k_train_test = obj.covar.calculateKernel(obj.x_data, test_Values);
            k_test_test = obj.covar.calculateKernel(test_Values, test_Values);

            avg_at_values = k_train_test'*obj.inv_covar_input*obj.y_data;
            
            covar_at_values = k_test_test - (k_train_test'*obj.inv_covar_input*k_train_test);

            variance = diag(covar_at_values);
            
            obj.memory = [avg_at_values; covar_at_values; variance];
            avg = avg_at_values;
        end
    end
end