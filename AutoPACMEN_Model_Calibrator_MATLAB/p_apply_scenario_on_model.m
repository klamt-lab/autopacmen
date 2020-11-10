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

function [model] = p_apply_scenario_on_model(model, scenario)
    % Applies the given scenario entry in JSON form on the given CNA model

    % Get the setup (i.e., the conditions for the scenario)
    setup_reaction_keys = fieldnames(scenario.setup);
    for j = 1:numel(setup_reaction_keys)
        % Get setup data for scenario
        reaction_name = setup_reaction_keys{j};
        bounds = eval("scenario.setup."+reaction_name);
        reaction_id = strmatch(reaction_name, model.reacID, 'exact');

        % Set lower and upper bounds \o/
        change_keys = fieldnames(bounds);
        for change_key_index = 1:numel(change_keys)
            change_key = change_keys{change_key_index};
            if change_key == "lower_bound"
                model.reacMin(reaction_id) = bounds.lower_bound;
            end
            if change_key == "upper_bound"
                model.reacMax(reaction_id) = bounds.upper_bound;
            end
        end
    end
end
