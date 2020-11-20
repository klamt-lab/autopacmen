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

function [best_prot_pool, start_prot_pool] = smoment_prot_pool_optimization(cna_model, scenarios, scenarios_matrix, lowest_value, highest_value)
    % AutoPACMEN protein pool model calibrator.
    %
    % Calibrates the protein pool so that it fits best with the given
    % scenarios. For a more detailed explanation on scenarios and this
    % function's context within AutoPACMEN, see AutoPACMEN's manual.
    %
    % Usage example:
    % [best_prot_pool, start_prot_pool] = smoment_prot_pool_optimization(cna_model,...
    %                                      scenarios,...
    %                                      scenarios_matrix,...
    %                                      max_change_factor)
    %
    % Arguments:
    % cna_model ~ A CellNetAnalyzer model, e.g. loaded with
    %             CNAsbmlModel2MFNetwork
    % scenarios ~ A struct from a scenarios JSON loaded with jsondecode.
    %             See the manual for more on this.
    % scenarios_matrix ~ The scenario dependency matrix. See the manual for
    %                    more on this.
    % lowest_value ~ The lowest possible intended protein pool.
    % highest_value ~ The highest possible intended protein pool.
    %
    % Outputs:
    % best_prot_pool ~ The final optimized prot_pool.
    % start_prot_pool ~ The unoptimized prot_pool. Useful for comparison
    %                   with start_prot_pool.
    
    % Set global variables used in fmincon
    global g_cna_model;
    g_cna_model = cna_model;
    global g_scenarios;
    g_scenarios = scenarios;
    global g_scenarios_matrix;
    g_scenarios_matrix = scenarios_matrix;
    
    % Get start prot pool value
    metabolite_index = strmatch('M_prot_pool', cna_model.specID, 'exact');
    reaction_index = strmatch('R_ER_pool_TG_', cna_model.specID, 'exact');
    start_prot_pool = cna_model.stoichMat(metabolite_index, reaction_index);
    
    % Set lower and upper bounds for fmincon
    lower_bounds = [lowest_value];
    upper_bounds = [highest_value];
    
    % Run optimization :D
    best_prot_pool = fmincon(@p_get_objective_value_for_prot_pool, start_prot_pool, [], [], [], [], lower_bounds, upper_bounds);
    
    % Delete global variables
    clear g_cna_model;
    clear g_reactions_to_change;
    clear g_scenarios;
    clear g_scenarios_matrix;
end
