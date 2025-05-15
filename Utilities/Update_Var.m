%% Update Variables with logical variable


function [Var_new] = Update_Var(Var, Exclusion)

    % Find all names and loop through them
    Variable_names = fieldnames(Var);
    
    % Preallocate space
    Var_new = struct;
    
    for x = 1:size(Variable_names)
        tmp_variable = Var.(Variable_names{x});
        if any(size(Exclusion, 1) == size(tmp_variable))
            Dim = find(size(Exclusion, 1) == size(tmp_variable));
            if Dim == 1
                Var_new.(Variable_names{x}) = tmp_variable(~Exclusion, :);
            elseif Dim == 2
                Var_new.(Variable_names{x}) = tmp_variable(:, ~Exclusion);
            elseif all(ismember(Dim, [1 2]))
                Var_new.(Variable_names{x}) = tmp_variable(~Exclusion, ~Exclusion);
            else
                disp('Might be a Problem here, check ! Code Update_Var, line 23 !')
                continue            
            end
        else
            Var_new.(Variable_names{x}) = tmp_variable;
        end
    end

end
