function [mcs]  = p_mcs_analysis(cnaModel,...
                                 outputFilename,...
                                 growthRateReactionName,...
                                 minGrowthRate,...
                                 substrateRateReactionName,...
                                 productRateReactionName,...
                                 minYield,...
                                 maxMCSsize)
    % Allow secretion
    cnaModel.reacMax(strmatch(productRateReactionName, cnaModel.reacID, 'exact')) = 1000;

    % Undesired fluxes
    T = zeros([2, cnaModel.numr]);
    % mu_min <= minGrowthRate h⁻1
    T(1, strmatch(growthRateReactionName, cnaModel.reacID)) = 1;
    t(1) = minGrowthRate;

    T(2, strmatch(substrateRateReactionName, cnaModel.reacID)) = minYield; % Substrate uptake has a negative flux :O
    T(2, strmatch(productRateReactionName, cnaModel.reacID)) = 1;
    t(2) = 0;

    t = t';

    D = [];
    d = [];
    % Desired fluxes
    D = zeros([2, cnaModel.numr]);
    % Product/Substrate >= minYield <=> minYield*Substrate - Product <= 0
    D(1, strmatch(substrateRateReactionName, cnaModel.reacID)) = -minYield; % Substrate uptake has a negative flux :O
    D(1, strmatch(productRateReactionName, cnaModel.reacID)) = -1;
    d(1) = 0;

    % -mu_min <= -minGrowthRate h⁻1
    D(2, strmatch(growthRateReactionName, cnaModel.reacID)) = -1;
    d(2) = -minGrowthRate;

    d = d';

    notknockable = [];
    maxMCSnum = 100000000;
    mcs = CNAregMCSEnumerator(cnaModel, T, t, D, d, notknockable, maxMCSnum, maxMCSsize, outputFilename);
end
