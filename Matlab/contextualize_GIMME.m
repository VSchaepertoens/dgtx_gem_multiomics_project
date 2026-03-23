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
exprThreshold = 0;  % since your data already contains zeros


%% 6.  Load experimental growth rates ---
growthTable = readtable(fullfile(PROJECT_ROOT,'data','rates','growth_rate.csv'));
% Optional: convert Condition & Window to strings if needed
growthTable.Condition = string(growthTable.Condition);
growthTable.Window    = string(growthTable.Window);

%% 7. Loop over models & condition-window columns

constrainedFolder = fullfile(PROJECT_ROOT, 'models/exp_constrained');

for m = 1:numel(modelNames)

    currentModelName = modelNames{m};
    fprintf('\n=== Processing model: %s ===\n', currentModelName);

    % Filter expression table for this model
    exprModel = exprTable(strcmp(exprTable.model, currentModelName), :);

    % Identify context columns
    contextCols = exprModel.Properties.VariableNames(3:end);

    % Convert context name to match constrained model naming
    %contextFile = strrep(contextCols, '_W', '_');

    % Biomass reaction
    objectiveRxn = biomassReactions(currentModelName);

   for i = 1:numel(contextCols)

    % Name in expression table
    contextExpr = contextCols{i};

    % Name used in model files
    contextFile = strrep(contextExpr, '_W', '_');

    fprintf('Processing context: %s\n', contextExpr);

    % --- Load correct constrained model ---
    constrainedModelFile = fullfile( ...
        constrainedFolder, ...
        sprintf('%s_%s_constrained.mat', currentModelName, contextExpr) ...
    );

    if ~isfile(constrainedModelFile)
        warning('Model not found: %s', constrainedModelFile);
        continue
    end

    model = readCbModel(constrainedModelFile);

    % % --- Inside the context loop, after loading the constrained model ---
    % % Get the experimental growth rate for this context
    % condStr = contextExpr; % e.g., 'TempShift_W1'
    % 
    % % Split contextExpr into condition and window
    % tokens = split(condStr,'_');  % {'TempShift','1'}
    % conditionName = tokens{1};
    % windowNum     = tokens{2};
    % 
    % % Find the matching row in growthTable
    % idx = strcmp(growthTable.Condition, conditionName) & strcmp(growthTable.Window, windowNum);
    % 
    % if any(idx)
    %     % Calculate lower bound
    %     expLB = growthTable.avg_growth__(idx) - growthTable.sd_growth__(idx);
    %     expLB = max(expLB, 0);  % avoid negative lower bound
    % 
    %     % Set biomass lower bound
    %     rxnID = findRxnIDs(model, objectiveRxn);
    %     model.lb(rxnID) = expLB;
    % 
    %     % Display current bounds
    %     fprintf('Setting biomass %s LB = %.5f, UB = %.1f\n', ...
    %             objectiveRxn, model.lb(rxnID), model.ub(rxnID));
    % else
    %     warning('No experimental growth rate found for %s window %s', conditionName, windowNum);
    % end

    % --- Display current biomass bounds ---
    biomassIdx = find(strcmp(model.rxns, objectiveRxn));
    if ~isempty(biomassIdx)
        fprintf('Biomass bounds for %s: LB = %.4f, UB = %.4f\n', ...
            objectiveRxn, model.lb(biomassIdx), model.ub(biomassIdx));
    else
        warning('Biomass reaction %s not found in model!', objectiveRxn);
    end

    % % --- Set threshold for GIMME (90% of experimental growth) ---
    % exprThreshold = 0.9 * expLB;  % 90% of experimental growth


    % --- Extract gene expression ---
    genes  = exprModel.EntrezID;
    values = exprModel.(contextExpr);   % original column name

    % --- Ensure proper types ---
    genes = cellstr(string(genes));   % convert to cell array of chars
    values = double(values);           % numeric vector

    % --- Create struct for COBRA ---
    geneExpStruct = struct();          % <-- make sure this is a struct
    geneExpStruct.gene  = genes;
    geneExpStruct.value = values;

    % --- Map expression to reactions ---
    rxnScores = mapExpressionToReactions(model, geneExpStruct);

    % Assign -1 to reactions not mapped to genes
    rxnScores(isnan(rxnScores)) = -1;

    % --- Run GIMME ---
    contextModel = GIMME(model, rxnScores, exprThreshold, objectiveRxn);
    % --- Quick stats ---
    fprintf('Reactions: original %d, context-specific %d, removed %d\n', ...
        numel(model.rxns), numel(contextModel.rxns), numel(model.rxns)-numel(contextModel.rxns));

    % --- Set model ID ---
    contextModel.id = sprintf('%s_%s_GIMME', currentModelName, contextExpr);

    % % --- Save ---
    % saveFile = fullfile(outFolder, [contextModel.id '.mat']);
    % writeCbModel(contextModel, 'fileName', saveFile, 'format', 'mat');
    % 
    % saveXmlFile = fullfile(outFolder, [contextModel.id '.xml']);
    % writeCbModel(contextModel, 'fileName', saveXmlFile, 'format', 'sbml');

     end
end

fprintf('\nAll models and contexts processed.\n');
