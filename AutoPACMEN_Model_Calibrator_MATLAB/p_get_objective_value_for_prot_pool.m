% Copyright 2019 PSB
%
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
%
%     http://www.apache.org/licenses/LICENSE-2.0
%
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.


function [f] = p_get_objective_value_for_prot_pool(x)
    % Returns f, the objective value (which shall be minimized) with the given
    % current state of optimized protein pool.

    % Get given global variables set in smoment_prot_pool_optimization.m
    global g_cna_model;
    global g_scenarios;
    global g_scenarios_matrix;

    objective_value = 0;

    % Get changed model
    changed_model = g_cna_model;

    % Change prot pool
    metabolite_index = strmatch('M_prot_pool', g_cna_model.specID, 'exact');
    reaction_index = strmatch('R_ER_pool_TG_', g_cna_model.specID, 'exact');
    g_cna_model.stoichMat(metabolite_index, reaction_index) = x;

    % Get values for each scenario line
    number_scenario_lines = size(g_scenarios_matrix, 1);
    fba_error = false;
    sum_scenario_weights = 0;
    for i = 1:number_scenario_lines
        line = g_scenarios_matrix(i,:);
        % Get scenario names \o/
        maximization_scenario_name = line(1);
        maximization_scenario_name = maximization_scenario_name{1};
        minimization_scenario_name = line(2);
        minimization_scenario_name = minimization_scenario_name{1};
        % Set maximization variables (used if subsequent minimization is
        % used :-)
        maximization_target_id = -1;
        maximization_result = NaN;
        % Run maximization scenario :D
        if maximization_scenario_name ~= "NA"
            % Apply scenario on model
            scenario = g_scenarios.(maximization_scenario_name);
            scenario_model = p_apply_scenario_on_model(changed_model, scenario);

            % Set model objective :-)
            maximization_target = scenario.target.reaction;
            maximization_target_id = strmatch(maximization_target, scenario_model.reacID, 'exact');
            scenario_model.objFunc(:) = 0;
            scenario_model.objFunc(maximization_target_id) = -1;

            % Perform FBA :D
            [~, success, ~, maximization_result] = CNAoptimizeFlux(scenario_model, [], [], 0, -1);
            maximization_result = maximization_result * -1;
            if (success == 0) || (isnan(maximization_result))
                disp("FBA error (maximization) D:");
                fba_error = true;
            end

            % Calculate objective value
            target_value = scenario.target.value;
            % objective_value_addition = abs( abs(maximization_result) - abs(target_value) );
            objective_value_addition = abs( maximization_result / target_value );
            if (objective_value_addition < 1.0)
                objective_value_addition = 1 / objective_value_addition;
            end
            objective_value_addition = objective_value_addition - 1;
            if strmatch('weight', fieldnames(scenario), 'exact')
                objective_value_addition = objective_value_addition * scenario.weight;
                sum_scenario_weights = sum_scenario_weights + scenario.weight;
            else
                sum_scenario_weights = sum_scenario_weights + 1;
            end
            objective_value = objective_value + abs(objective_value_addition);
        end

        % Run minimization scenario :D
        if minimization_scenario_name ~= "NA"
            scenario = g_scenarios.(minimization_scenario_name);
            scenario_model = p_apply_scenario_on_model(changed_model, scenario);
            if maximization_target_id ~= -1
                scenario_model.reacMin(maximization_target_id) = maximization_result * 0.9999;
                scenario_model.reacMax(maximization_target_id) = maximization_result * 1.0001;
            end

            % Set model objective :-)
            minimization_target = scenario.target.reaction;
            minimization_target_id = strmatch(minimization_target, scenario_model.reacID, 'exact');
            scenario_model.objFunc(:) = 0;
            scenario_model.objFunc(minimization_target_id) = -1;

            % Perform FBA :D
            [~, success, ~, minimization_result] = CNAoptimizeFlux(scenario_model, [], [], 0, -1);
            if (success == 0) || (isnan(minimization_result))
                disp("FBA error (minimization) D:");
                fba_error = true;
            end

            % Calculate objective value
            target_value = scenario.target.value;
            % objective_value_addition = abs( abs(minimization_result) - abs(target_value) );
            objective_value_addition = abs( minimization_result / target_value );
            if (objective_value_addition < 1.0)
                objective_value_addition = 1 / objective_value_addition;
            end
            objective_value_addition = objective_value_addition - 1;
            if strmatch('weight', fieldnames(scenario), 'exact')
                objective_value_addition = objective_value_addition * scenario.weight;
                sum_scenario_weights = sum_scenario_weights + scenario.weight;
            else
                sum_scenario_weights = sum_scenario_weights + 1;
            end
            objective_value = objective_value + abs(objective_value_addition);
        end
    end

    % Return objective value :D
    objective_value = objective_value / sum_scenario_weights;
    if ~fba_error
        f = objective_value;
        disp(objective_value)
    else
        f = 100000;
    end
    % sprintf('%.12f', objective_value)
end