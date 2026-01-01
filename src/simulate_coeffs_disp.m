function simulate_coeffs_disp(coeffs, struct_name)
    if nargin < 2
        struct_name  = "coeffs";
    end

    keys = fieldnames(coeffs);
    for i = 1:length(keys)
        key = keys{i};
        if length(coeffs.(key)) > 1 
            if size(coeffs.(key), 1) == 1
                disp(struct_name + "." + key + " = [" + sprintf('%0.4f ', coeffs.(key)) + "];")
            else
                disp(struct_name + "." + key + " = [" + sprintf('%0.4f ', coeffs.(key)) + "]';")
            end
        else
            disp(struct_name  + "." + key + " = " + sprintf('%0.4f ', coeffs.(key)) + ";")
        end
    end
    disp(" ")
end
