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

function [best_kcats, start_kcats] = smoment_kcat_optimization(cna_model, reactions_to_change, scenarios, scenarios_matrix, max_change_factor)
    % AutoPACMEN kcat model calibrator.
    %
    % Calibrates the protein pool so that it fits best with the given
    % scenarios. For a more detailed explanation on scenarios and this
    % function's context within AutoPACMEN, see AutoPACMEN's manual.
    %
    % Usage example:
    % [best_kcats, start_kcats] = smoment_kcat_optimization(cna_model,...
    %                                                       reactions_to_change,...
    %                                                       scenarios,...
    %                                                       scenarios_matrix,...
    %                                                       max_change_factor)
    %
    % Arguments:
    % cna_model ~ A CellNetAnalyzer model, e.g. loaded with
    %             CNAsbmlModel2MFNetwork
    % reactions_to_change ~ A list of reactions for which their kcats
    %                       shall be optimized.
    % scenarios ~ A struct from a scenarios JSON loaded with jsondecode.
    %             See the manual for more on this.
    % scenarios_matrix ~ The scenario dependency matrix. See the manual for
    %                    more on this.
    % max_change_factor ~ Number; The kcats will be restriced to be
    %                     minimally max_change_factor times lower or
    %                     max_change_factor times higher than their start
    %                     value.
    %
    % Outputs:
    % best_kcats ~ A list of optimized MW/kcat values in the order of
    %              the reactions given in reactions_to_change.
    % start_kcats ~ The unoptimized MW/kcat values in the order of
    %               the reactions given in reactions_to_change. Useful in
    %               order to compare them with best_kcats.
    
    % Set global variables used in fmincon
    global g_cna_model;
    g_cna_model = cna_model;
    global g_reactions_to_change
    g_reactions_to_change = reactions_to_change;
    global g_scenarios;
    g_scenarios = scenarios;
    global g_scenarios_matrix;
    g_scenarios_matrix = scenarios_matrix;
    
    % Get start kcats
    start_kcats = [];
    metabolite_index = strmatch('M_prot_pool', cna_model.specID, 'exact');
    for i = 1:numel(reactions_to_change)
        reaction_id = strmatch(reactions_to_change{i}, cna_model.reacID, 'exact');
        start_kcats(i) = cna_model.stoichMat(metabolite_index, reaction_id);
    end
    
    % Set lower and upper bounds for fmincon
    lower_bounds = [];
    upper_bounds = [];
    for i = 1:numel(start_kcats)
        lower_bounds(i) = (start_kcats(i)) * max_change_factor;
        upper_bounds(i) = (start_kcats(i)) / max_change_factor;
    end
    
    % Run optimization :D
    best_kcats = fmincon(@p_get_objective_value_for_kcats, start_kcats, [], [], [], [], lower_bounds, upper_bounds);
    
    % Delete global variables
    clear g_cna_model;
    clear g_reactions_to_change;
    clear g_scenarios;
    clear g_scenarios_matrix;
end
