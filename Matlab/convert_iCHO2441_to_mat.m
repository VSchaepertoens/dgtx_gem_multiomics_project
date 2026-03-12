%% 1. Setup
% Add COBRA Toolbox paths if not already added
initCobraToolbox(false)   % false = do not update
%%
% Load XML/SBML model
modelXml = 'models/original/iCHO2441.xml';
model = readCbModel(modelXml);

% Save as MAT
saveFile = 'models/original/iCHO2441.mat';
writeCbModel(model, 'fileName', saveFile, 'format', 'mat');

fprintf('Model saved as MAT: %s\n', saveFile);