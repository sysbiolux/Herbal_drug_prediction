WORKING_DIR = ''
cd(WORKING_DIR)
cd('./Herbal_drug_prediction/scripts/')
addpath(genpath('./Herbal_drug_prediction/scripts'))
% Calculate the rxn presence rate for each pathway for each model

%%loading generic model
genericModel = readCbModel('../inputs/Recon3DModel_301.mat');
genericModel.subSystems = string(genericModel.subSystems);
%creating an empty table
final_table = [];
%read all .mat files inside folder
files = dir('../models/models after removal of unused genes/*.mat');
files = {files.name}';
for i=1:length(files)
    file = strcat('../models/models after removal of unused genes/',files(i));
    contextModel = readCbModel(file{1});
    %read the model
    % call the function rxn_per_pathway
    contextModel.subSystems = string(contextModel.subSystems);
    [rxn_per_pathway, ~] = contextModelEnrichment(contextModel,genericModel);
    drug_name= cellstr(repmat(files(i), numel(rxn_per_pathway.pathway),1));
    rxn_per_pathway_name = [drug_name rxn_per_pathway] ;
    final_table = [final_table;rxn_per_pathway_name];
end

final_table.Properties.VariableNames = ({'drug', 'pathway', 'Rxn_per_pathway'});
strg_to_replace = {'Four mixture(Tanshinone IIA:Salvianic acid A sodium:Protocatechuic aldehyde:Salvianolic acid B=1:1:1:1)', '(+)2-(1-hydroxyl-4-oxocyclohexyl) ethyl caffeate', '_model.mat'};
new_names = {'Four_mixture', 'Ethyl_caffeate', ''};
final_table.drug = replace(final_table.drug, strg_to_replace, new_names) ;
Models_rxn_per_pathway = final_table;
writetable(Models_rxn_per_pathway,'../results/Models_rxn_per_pathway.csv');

%% count the total number of reactions for each pathways in Recon3

pathway_count = num2cell(histcounts(categorical(string(genericModel.subSystems))));
pathway_names = (unique(string(genericModel.subSystems)));
count_pathway = [pathway_names,pathway_count'];
count_pathway= cellstr(count_pathway);
count_pathway = cell2table(count_pathway);
count_pathway.Properties.VariableNames = [{'pathway','count'}];
writetable(count_pathway, '../results/Recon3d_pathway_total_counts.csv')

