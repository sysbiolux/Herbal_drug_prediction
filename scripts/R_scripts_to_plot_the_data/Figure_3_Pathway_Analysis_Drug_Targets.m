WORKING_DIR = ''
cd(WORKING_DIR)
cd('./Herbal_drug_prediction/scripts/')
addpath(genpath('./Herbal_drug_prediction/scripts'))
%% Calculate the drugs target statistics for the herbal and cancer drugs
if exist('../inputs/dico_201911.mat.gz')>=1 % Extracting the compressed dictionary
    gunzip('../inputs/dico_201911.mat.gz')
end
load('../inputs/dico_201911.mat');
% How many unique targets
% How many unique rxns for the drug targets
% How many unique pathways for the drug targets
target_db = load('../inputs/Drug_target_all_databases_herbs.mat');
target_db = target_db.drug_target_all_database_herbs;
recon3d = readCbModel('../inputs/Recon3DModel_301.mat');
genes = recon3d.genes;
for i=1:numel(genes)
    x = strsplit(genes{i,1},'.');
    genes{i} = x{1};
end
recon3d.genes = genes;

% How many unique targets = 474
numel(unique(target_db.entrez_id))

% How many unique metabolic targets = 96
target_genes = intersect(string(unique(target_db.entrez_id)),recon3d.genes);
numel(unique(target_genes))

%% How many herbal drugs across the three databases =44
numel(unique(target_db.name))


%% How many herbal drugs with metabolic target =35
numel(unique(target_db(ismember(string(target_db.entrez_id),target_genes),4)))

% How many unique rxns for the drug targets = 890
rxns_associated_genes = recon3d.rxns(~strcmp(recon3d.rules,''));
rxns_associated_genes_idxs = find(~strcmp(recon3d.rules,''));

target_genes_idx =find(ismember(recon3d.genes,target_genes));
recon3d.genes(target_genes_idx)
[r_genes,~]=find(recon3d.rxnGeneMat(:, target_genes_idx));
numel(unique(intersect(r_genes,rxns_associated_genes_idxs)))

% How many unique pathways for the drug targets =49
genes_ss=recon3d.subSystems(unique(r_genes));
numel(unique(cell2table( genes_ss)))

%Mapping gene entrez to symbols
target_db.symbol = string([1:size(target_db.entrez_id,1)]');
target_db.is_metabolic = ismember(string(target_db.entrez_id),recon3d.genes);

for i=1:size(target_db,1)
    entrez= target_db{i,5};
    symbol = dico(find(ismember(dico.ENTREZ,string(entrez))),5);
    target_db(i,6) = unique(table2cell(symbol));
end

% Table of druggable pathways
herbal_pathways = unique(cell2table( genes_ss));
herbal_pathways.total_rxns = [1:size(herbal_pathways.genes_ss,1)]';
herbal_pathways.total_genes = [1:size(herbal_pathways.genes_ss,1)]';
herbal_pathways.targeted_rxns = [1:size(herbal_pathways.genes_ss,1)]';
herbal_pathways.targeted_genes = [1:size(herbal_pathways.genes_ss,1)]';
herbal_pathways.targeted_drugs = [1:size(herbal_pathways.genes_ss,1)]';
herbal_pathways.targeted_drugs_names = repmat("NA", size(herbal_pathways.genes_ss,1),1);
herbal_pathways.targeted_genes_names = repmat("NA", size(herbal_pathways.genes_ss,1),1);
pathways = cell2table(recon3d.subSystems);
for p=1:size(herbal_pathways.genes_ss,1)
    pathway = herbal_pathways.genes_ss{p};
    pathway
    pathway_idx = find(ismember(pathways.Var1, pathway));
    pathway_idx = intersect(rxns_associated_genes_idxs,pathway_idx);
    herbal_pathways.total_rxns(p) = numel(unique(recon3d.rxns(pathway_idx))); %numel(intersect(recon3d.rxns(pathway_idx),rxns_associated_genes));%;%
    pathway_rxns = recon3d.rxns(pathway_idx);
    [~,r_genes_p]=find(recon3d.rxnGeneMat(pathway_idx, :));
    herbal_pathways.total_genes(p) = numel(unique(r_genes_p));
    herbal_pathways.targeted_rxns(p) = numel(intersect(unique(r_genes),pathway_idx));
    targeted_genes_p_idx = intersect(unique(r_genes_p),target_genes_idx);
    herbal_pathways.targeted_genes(p) = numel(targeted_genes_p_idx);
    herbal_pathways.targeted_genes_names(p) = string(strjoin(unique(recon3d.genes(targeted_genes_p_idx)),'; '));

    shared_drugs_in_p = unique(target_db(ismember(string(target_db.entrez_id),string(recon3d.genes(targeted_genes_p_idx))),4));
    herbal_pathways.targeted_drugs(p) = numel(unique(shared_drugs_in_p));
    herbal_pathways.targeted_drugs_names(p) = string(strjoin(shared_drugs_in_p.name,'; '));
end
herbal_pathways.target_gene_ratio =herbal_pathways.targeted_genes./ herbal_pathways.total_genes;
herbal_pathways.target_rxn_ratio =herbal_pathways.targeted_rxns./ herbal_pathways.total_rxns;

pathways_of_interest= {'Androgen, estrogen synthesis and metabolism','Arachidonic acid metabolism','Drug metabolism','Steroid metabolism','Linoleate metabolism'};
pathways_idx = find(ismember(herbal_pathways.genes_ss,pathways_of_interest));
pathways_drugs = herbal_pathways.targeted_drugs_names(pathways_idx);
pathways_drugs = strjoin(pathways_drugs,'; ');
numel(unique(strsplit(pathways_drugs,'; ')))

writetable(herbal_pathways,'../results/Summary_Herbal_drugs_targets.csv');


%% For pharmaceutical drugs
target_db = load('../inputs/GeneTargetPharma.mat');
target_db = target_db.GeneTargetPharma;
% How many unique targets = 335
numel(unique(target_db.ENTREZ))

% How many unique metabolic targets = 82
target_genes = intersect(string(unique(target_db.ENTREZ)),recon3d.genes);
numel(target_genes)

%% How many pharma drugs with metabolic target =26
numel(unique(target_db(ismember(string(target_db.ENTREZ),target_genes),2)))

% How many unique rxns for the drug targets = 1207
target_genes_idx =find(ismember(recon3d.genes,target_genes));
recon3d.genes(target_genes_idx)
[r_genes,~]=find(recon3d.rxnGeneMat(:, target_genes_idx));
numel(unique(r_genes))

% How many unique pathways for the drug targets =45
genes_ss=recon3d.subSystems(unique(r_genes));
numel(unique(cell2table( genes_ss)))

%Mapping gene entrez to symbols
target_db.symbol = string([1:size(target_db.ENTREZ,1)]');
target_db.is_metabolic = ismember(string(target_db.ENTREZ),recon3d.genes);

for i=1:size(target_db,1)
    entrez= target_db{i,2};
    symbol = dico(find(ismember(dico.ENTREZ,string(entrez))),5);
    target_db(i,3) = unique(table2cell(symbol));
end

% Table of druggable pathways
pharma_pathways = unique(cell2table( genes_ss));
pharma_pathways.total_rxns = [1:size(pharma_pathways.genes_ss,1)]';
pharma_pathways.total_genes = [1:size(pharma_pathways.genes_ss,1)]';
pharma_pathways.targeted_rxns = [1:size(pharma_pathways.genes_ss,1)]';
pharma_pathways.targeted_genes = [1:size(pharma_pathways.genes_ss,1)]';
pharma_pathways.targeted_drugs = [1:size(pharma_pathways.genes_ss,1)]';
pharma_pathways.targeted_drugs_names = repmat("NA", size(pharma_pathways.genes_ss,1),1);
pharma_pathways.targeted_genes_names = repmat("NA", size(pharma_pathways.genes_ss,1),1);

rxns_associated_genes = recon3d.rxns(~strcmp(recon3d.rules,''));
rxns_associated_genes_idxs = find(~strcmp(recon3d.rules,''));

pathways = cell2table(recon3d.subSystems);
for p=1:size(pharma_pathways.genes_ss,1)
    pathway = pharma_pathways.genes_ss{p};
    %pathway
    pathway_idx = find(ismember(pathways.Var1, pathway));
    pathway_idx = intersect(rxns_associated_genes_idxs,pathway_idx);
    pharma_pathways.total_rxns(p) = numel(recon3d.rxns(pathway_idx));
    pathway_rxns = recon3d.rxns(pathway_idx);
    [~,r_genes_p]=find(recon3d.rxnGeneMat(pathway_idx, :));
    pharma_pathways.total_genes(p) = numel(unique(r_genes_p));
    pharma_pathways.targeted_rxns(p) = numel(intersect(unique(r_genes),pathway_idx));
    targeted_genes_p_idx = intersect(unique(r_genes_p),target_genes_idx);
    pharma_pathways.targeted_genes(p) = numel(targeted_genes_p_idx);

    shared_drugs_in_p = unique(target_db(ismember(string(target_db.ENTREZ),string(recon3d.genes(targeted_genes_p_idx))),1));
    pharma_pathways.targeted_drugs(p) = numel(shared_drugs_in_p);
    pharma_pathways.targeted_drugs_names(p) = string(strjoin(shared_drugs_in_p.DrugName,'; '));
    pharma_pathways.targeted_genes_names(p) = string(strjoin(unique(recon3d.genes(targeted_genes_p_idx)),'; '));
end
pharma_pathways.target_gene_ratio =pharma_pathways.targeted_genes./ pharma_pathways.total_genes;
pharma_pathways.target_rxn_ratio =pharma_pathways.targeted_rxns./ pharma_pathways.total_rxns;

writetable(pharma_pathways,'../results/Summary_Pharmaceutical_drugs_targets.csv');
