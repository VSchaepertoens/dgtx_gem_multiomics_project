%% GIMME Contextualization for CHO GEMs
% This script runs GIMME on CHO GEMs for each condition-window using
% preprocessed transcriptomics data.

%% 1. Setup
% Add COBRA Toolbox paths if not already added
initCobraToolbox(false)   % false = do not update

%% 2. Project paths
% Project root
PROJECT_ROOT = pwd; 
outFolder = fullfile(PROJECT_ROOT, 'models/contextualized/GIMME');
if ~exist(outFolder, 'dir')
    mkdir(outFolder);
end

%% 3. Define original GEMs
modelFiles = { ...
    'models/original/iCHO3K_cho_prod_generic_unblocked.mat', ...
    'models/original/iCHOv1.mat' ...
    'models/original/iCHO2441.mat', ...
    };

modelNames = { ...
    'iCHO3K_cho_prod_generic_unblocked', ...
    'iCHOv1' ...
    'iCHO2441', ...
    };

biomassReactions = containers.Map( ...
    modelNames, ...
    {'biomass_cho_prod', 'BIOMASS_cho_producing', 'biomass_cho_prod'} ...
        );


%% 4. Load transcriptomics data
% CSV format: EntrezID, Constant_W1, Constant_W2, ..., TempShift_W5
exprTable = readtable('data/transcriptomics/charrun_FBA_RNA3_all_V02_20260311_VS.csv');

% Optional: convert EntrezID column to strings to match COBRA Toolbox gene IDs
exprTable.EntrezID = string(exprTable.EntrezID);

% --- Rename context columns to match constrained model filenames ---
varNames = exprTable.Properties.VariableNames;

for i = 1:numel(varNames)
    name = varNames{i};

    % Skip non-context columns
    if strcmp(name, 'model') || strcmp(name, 'EntrezID')
        continue
    end

    newName = name;

    % Replace '_W' → '_'
    newName = strrep(newName, '_W', '_');

    % Replace 'TempShift' → 'Temp shifted'
    newName = strrep(newName, 'TempShift', 'Temp. shifted');

    % Apply new name
    exprTable.Properties.VariableNames{i} = newName;
end

%% 5. Define GIMME parameters
% Threshold for keeping reactions (can be >0 or another value)
exprThreshold = 1;  


% %% 6.  Load experimental growth rates ---
% growthTable = readtable(fullfile(PROJECT_ROOT,'data','rates','growth_rate.csv'));
% % Optional: convert Condition & Window to strings if needed
% growthTable.Condition = string(growthTable.Condition);
% growthTable.Window    = string(growthTable.Window);

%% 
function [contextModel, diagStruct] = processGIMMEContext(model, exprModel, contextExpr, exprThreshold, biomassFunc, modelName)
% PROCESSGIMMECONTEXT Extracts context-specific model using GIMME
%
% Inputs:
%   model       - COBRA model structure
%   exprModel   - table with gene expression for current model
%   contextExpr - name of the context column in exprModel
%   exprThreshold - threshold for GIMME (RPKM/FPKM > 1)
%   biomassFunc   - function handle or map to get correct biomass reaction
%
% Output:
%   contextModel - GIMME-constrained context-specific model
%   diagStruct - Diagnostics for each context

    % --- Extract gene expression ---
    genes  = exprModel.EntrezID;
    values = exprModel.(contextExpr);
    
    % --- Ensure proper types ---
    genes = cellstr(string(genes));   
    values = double(values);          
    
    % --- Create struct for COBRA ---
    geneExpStruct = struct();         
    geneExpStruct.gene  = genes;
    geneExpStruct.value = values;

    % --- Map expression to reactions ---
    [rxnScores, parsedGPR, gene_used] = mapExpressionToReactions(model, geneExpStruct);
    rxnScores(isnan(rxnScores)) = -1;

    % --- Diagnostics ---
    numRxns = numel(model.rxns);
    numMapped = sum(rxnScores >= 0);
    numUnmapped = numRxns - numMapped;
    fracAboveThreshold = sum(rxnScores >= exprThreshold) / numRxns;
    overlapGenes = intersect(model.genes, geneExpStruct.gene);

    diagStruct = struct();
    diagStruct.modelName = modelName;
    diagStruct.context = contextExpr;
    diagStruct.numRxns = numRxns;
    diagStruct.numMapped = numMapped;
    diagStruct.numUnmapped = numUnmapped;
    diagStruct.fracAboveThreshold = fracAboveThreshold;
    diagStruct.minRxnScore = min(rxnScores);
    diagStruct.medianRxnScore = median(rxnScores);
    diagStruct.maxRxnScore = max(rxnScores);
    diagStruct.numGenesWithExpr = numel(overlapGenes);
    diagStruct.numGenesTotal = numel(model.genes);

    fprintf('Reactions mapped: %d / %d, unmapped: %d\n', numMapped, numRxns, numUnmapped);
    fprintf('Reaction scores: min=%.4f, median=%.4f, max=%.4f\n', ...
        diagStruct.minRxnScore, diagStruct.medianRxnScore, diagStruct.maxRxnScore);
    fprintf('Fraction of reactions above threshold (%.2f): %.2f%%\n', ...
        exprThreshold, fracAboveThreshold*100);
    fprintf('Genes with expression data: %d / %d\n', ...
        diagStruct.numGenesWithExpr, diagStruct.numGenesTotal);

    % --- Set the correct objective ---
    objectiveRxn = biomassFunc();
    model = changeObjective(model, objectiveRxn);
    fprintf('Objective set to %s\n', objectiveRxn);

    % --- Run GIMME ---
    contextModel = GIMME(model, rxnScores, exprThreshold);
    fprintf('Reactions: original %d, context-specific %d, removed %d\n', ...
        numRxns, numel(contextModel.rxns), numRxns - numel(contextModel.rxns));

    % --- Set model ID ---
    contextModel.id = sprintf('%s_%s_GIMME', modelName, contextExpr);

    % --- Add GIMME stats to diagnostics ---
    diagStruct.numRxnsContext = numel(contextModel.rxns);
    diagStruct.numRxnsRemoved = numRxns - numel(contextModel.rxns);

end

%% 7. Loop over models & condition-window columns

constrainedFolder = fullfile(PROJECT_ROOT, 'models/exp_constrained');

allDiags = [];  % initialize struct array

for m = 1:numel(modelNames)
    currentModelName = modelNames{m};
    fprintf('\n=== Processing model: %s ===\n', currentModelName);

    exprModel = exprTable(strcmp(exprTable.model, currentModelName), :);
    contextCols = exprModel.Properties.VariableNames(3:end);

    for i = 1:numel(contextCols)
        contextExpr = contextCols{i};
        fprintf('Processing context: %s\n', contextExpr);

        constrainedModelFile = fullfile(constrainedFolder, ...
            sprintf('%s_%s_constrained.mat', currentModelName, contextExpr));

        if ~isfile(constrainedModelFile)
            warning('Model not found: %s', constrainedModelFile);
            continue
        end

        model = readCbModel(constrainedModelFile);

        % Call the function and get diagnostics
        [contextModel, diagStruct] = processGIMMEContext(model, exprModel, contextExpr, exprThreshold, ...
            @() biomassReactions(currentModelName), currentModelName);

        %    --- Save ---
        saveFile = fullfile(outFolder, [contextModel.id '.mat']);
        writeCbModel(contextModel, 'fileName', saveFile, 'format', 'mat');
        saveXmlFile = fullfile(outFolder, [contextModel.id '.xml']);
        writeCbModel(contextModel, 'fileName', saveXmlFile, 'format', 'sbml');

        % Append diagnostics to array
        allDiags = [allDiags; diagStruct];  
    end
end

% --- Convert to table and save ---
diagTable = struct2table(allDiags);
writetable(diagTable, fullfile(outFolder, 'GIMME_diagnostics.csv'));
fprintf('\nAll diagnostics saved to GIMME_diagnostics.csv\n');