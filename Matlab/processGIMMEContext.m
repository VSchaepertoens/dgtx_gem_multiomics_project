function [contextModel, diagStruct, rxnScores,  parsedGPR, gene_used] = processGIMMEContext(model, exprModel, contextExpr, exprThreshold, biomassFunc, modelName)
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
    [rxnScores, parsedGPR, gene_used] = mapExpressionToReactions(model, geneExpStruct, true);
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