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

function [best_kcats, start_kcats] = kcat_optimization_run_fmincon_ec_model_2019_06_25()
    % This is the optimization run of iJO1366*C (iJO1366* with manually changes protein pool
    % and manual changes for anaerobic conditions) for the scenarios given in "optimization_scenarios.json".

    % Load Orth model with kcats :D
    [cna_model, ~] = CNAsbmlModel2MFNetwork('/.../iJO1366_sMOMENT_2019_06_25_STANDARD_EXCHANGE_SCENARIO_MANUAL_CHANGES.xml');
    % Fitted for 25 growth scenarios
    cna_model.reacMax(strmatch('R_ER_pool_TG_', cna_model.reacID, 'exact')) = .095;

    % Set arguments for optimization \o/
    % Reactions of which the kcats shall be changed :-)
    reactions_to_change = [
        % Acetate
        "R_FLDR2",
        "R_ACKr_TG_forward",
        "R_ACtex_TG_forward",
        "R_POR5_TG_reverse",
        "R_PTAr_TG_reverse",
        % Glycerol
        "R_GLYK",
        "R_TRDR",
        "R_G3PD2_TG_forward",
        "R_GLYCtex_TG_forward",
        "R_GTHOr_TG_forward",
        % Oxoglutarate
        "R_AKGt2rpp_TG_forward",
        "R_AKGtex_TG_forward",
        % L-Alanine
        "R_ASPT",
        "R_DAAD",
        "R_PROD2",
        "R_ALAR_TG_forward",
        "R_ALATA_L_TG_forward",
        "R_ALAtex_TG_forward",
        "R_GLUDy_TG_forward",
        "R_VALTA_TG_forward",
        % Pyruvate
        "R_PYRtex_TG_forward",
        % Fructose
        "R_FRUK",
        "R_FRUpts2pp",
        "R_FRUptspp",
        "R_FRUtex_TG_forward",
        % Guanosine
        "R_ALLTAMH",
        "R_ALLTN",
        "R_GMPR",
        "R_GSNK",
        "R_GSNt2pp",
        "R_GUAD",
        "R_UGLYCH",
        "R_XAND",
        "R_GSNtex_TG_forward",
        "R_PPM_TG_forward",
        "R_PRPPS_TG_reverse",
        "R_PUNP3_TG_forward",
        % N-acetylglucosamine
        "R_ACGAptspp",
        "R_AGDC",
        "R_ACGAtex_TG_forward",
        % Fumarate
        "R_FUMt2_2pp",
        "R_FUMt2_3pp",
        "R_FUMtex_TG_forward",
        % Ribose
        "R_PGCD",
        "R_PSP_L",
        "R_RBK",
        "R_RIBabcpp",
        "R_GHMT2r_TG_forward",
        "R_MTHFD_TG_forward",
        "R_RIBtex_TG_forward",
        % L-Lactate
        "R_L_LACt2rpp_TG_forward",
        "R_L_LACtex_TG_forward",
        % Gluconate
        "R_GLCNt2rpp_TG_forward",
        "R_GLCNtex_TG_forward",
        % Sorbitol
        "R_SBTptspp",
        "R_SBTPD_TG_forward",
        "R_SBTtex_TG_forward",
        % Glucosamine
        "R_GAMptspp",
        "R_GAMtex_TG_forward",
        % L-Malate
        "R_MALt2_2pp",
        "R_MALt2_3pp",
        "R_MALtex_TG_forward",
        % Succinate
        "R_SUCCt2_2pp",
        "R_SUCCt2_3pp",
        "R_SUCCtex_TG_forward",
        % Maltose
        "R_AMALT2",
        "R_AMALT3",
        "R_AMALT4",
        "R_MALTabcpp",
        "R_MALTtexi",
        "R_MLTP1_TG_forward",
        "R_MLTP2_TG_forward",
        "R_MLTP3_TG_forward",
        % Glucose 6-Phosphate
        "R_G6Pt6_2pp",
        "R_G6Ptex_TG_forward",
        "R_PItex_TG_reverse",
        % Mannitol
        "R_MNLptspp",
        "R_M1PD_TG_forward",
        "R_MNLtex_TG_forward",
        % Trehalose
        "R_TRE6PH",
        "R_TREHpp",
        "R_TREtex_TG_forward",
        % Xylose
        "R_XYLK",
        "R_XYLt2pp",
        "R_RPE_TG_reverse",
        "R_RPI_TG_reverse",
        "R_XYLI1_TG_forward",
        "R_XYLtex_TG_forward",
        % Mannose
        "R_MANptspp",
        "R_MAN6PI_TG_forward",
        "R_MANtex_TG_forward",
        % Galactose
        "R_GALt2pp",
        "R_GALKr_TG_forward",
        "R_GALtex_TG_forward",
        "R_UDPG4E_TG_reverse",
        "R_UGLT_TG_forward"
    ];

    % Read scenarios JSON
    text = fileread('./0FminconOptimization/optimization_scenarios.json');
    scenarios = jsondecode(text);
    % Get scenarios maximization/minimization matrix
    scenarios_matrix = {
        "Glucose_Aerobic_Growth", "NA";
        "Glucose_Anaerobic_Growth", "NA";
        "Acetate", "NA";
        "Glycerol", "NA";
        "Oxoglutarate", "NA";
        "L_Alanine", "NA";
        "Pyruvate", "NA";
        "Fructose", "NA";
        "Guanosine", "NA";
        "N_acetylglucosamine", "NA";
        "Fumarate", "NA";
        "Ribose", "NA";
        "L_Lactate", "NA";
        "Gluconate", "NA";
        "Sorbitol", "NA";
        "Glucosamine", "NA";
        "L_Malate", "NA";
        "Succinate", "NA";
        "Maltose", "NA";
        "Glucose_6_Phosphate", "NA";
        "Mannitol", "NA";
        "Trehalose", "NA";
        "Xylose", "NA";
        "Mannose", "NA";
        "Galactose", "NA";
    };

    max_change_factor = 50;
    [best_kcats, start_kcats] = moment_optimization_fmincon(cna_model, reactions_to_change, scenarios, scenarios_matrix, max_change_factor);
    disp(best_kcats)
    best_kcats = best_kcats';
    save('./0FminconOptimization/best_values_2019_06_25_manual_changes_change_factor_50.mat','best_kcats');
    start_kcats = start_kcats';
    save('./0FminconOptimization/start_values_2019_06_25_manual_changes_change_factor_50.mat','start_kcats');
end
