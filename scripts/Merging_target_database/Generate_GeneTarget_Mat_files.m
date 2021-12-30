%% Generate 2 tables as mat files (GeneTargetPharma & GeneTargetHerbal)
%% for drug deletion prediction that have metabolic target genes and entrez and drug names
WORKING_DIR = ''
cd(WORKING_DIR)
cd('./Herbal_drug_prediction/scripts/')

load('../inputs/dico_201911.mat');
writematrix('../inputs/dico_201911.csv','dico');

recon_model = readCbModel('../inputs/Recon3DModel_301.mat');
genes = recon_model.genes;
for i=1:numel(genes)
    x = strsplit(genes{i,1},'.');
    genes{i} = x{1};
end
recon_model.genes = genes;

% GeneTargetHerbal.mat
target_db = readtable('../inputs/drug_target_all_database_only_human_cancer_drugs.csv', 'Delimiter', ',');
numel(unique(target_db.name))
target_db_metabolic = target_db(find(ismember(string(target_db.entrez_id),string(recon_model.genes))),:);
target_db_metabolic = unique(target_db_metabolic(:,ismember(target_db_metabolic.Properties.VariableNames,["name","entrez_id"])));
target_db_metabolic.Properties.VariableNames = {'DrugName','ENTREZ'};
GeneTargetPharma = target_db_metabolic;
GeneTargetPharma.ENTREZ = string(GeneTargetPharma.ENTREZ);
save('../inputs/GeneTargetPharma.mat','GeneTargetPharma');

% Drug_target_all_databases_herbs.mat
target_db = readtable('../inputs/drug_target_all_database_only_human_herb_drugs.csv', 'Delimiter', ',');
target_db.x___3 = [];
drug_target_all_database_herbs = target_db;
save('../inputs/Drug_target_all_databases_herbs.mat','drug_target_all_database_herbs');

numel(unique(target_db.name))
target_db_metabolic = target_db(find(ismember(string(target_db.entrez_id),string(recon_model.genes))),:);
target_db_metabolic = unique(target_db_metabolic(:,ismember(target_db_metabolic.Properties.VariableNames,["name","entrez_id"])));
target_db_metabolic.Properties.VariableNames = {'DrugName','ENTREZ'};
GeneTargetHerbal = target_db_metabolic;
save('../inputs/GeneTargetHerbal.mat','GeneTargetHerbal');
