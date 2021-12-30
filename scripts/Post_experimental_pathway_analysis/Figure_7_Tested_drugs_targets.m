WORKING_DIR = ''
cd(WORKING_DIR)
cd('./Herbal_drug_prediction/scripts/')

%% create a table of the number of drugs and rxns for each gene-pathway
% For pharmaceutical drugs
target_db = load('../inputs/GeneTargetPharma.mat');
target_db = target_db.GeneTargetPharma;
recon_model = readCbModel('../inputs/Recon3DModel_301.mat');
load('../inputs/dico_201911.mat')
genes = recon_model.genes;
for i=1:numel(genes)
    x = strsplit(genes{i,1},'.');
    genes{i} = x{1};
end
recon_model.genes = genes;
target_genes = intersect(string(unique(target_db.ENTREZ)),recon_model.genes);
target_genes_idx =find(ismember(recon_model.genes,target_genes));
recon_model.genes(target_genes_idx)
% all rxns of the all cancer drug targets
[r_genes,~]=find(recon_model.rxnGeneMat(:, target_genes_idx));
numel(unique(r_genes))

% How many unique pathways for the drug targets =45
genes_ss=recon_model.subSystems(unique(r_genes));
genes_ss = unique(string(genes_ss));
numel(unique( genes_ss))

%Mapping gene entrez to symbols
target_db.symbol = string([1:size(target_db.ENTREZ,1)]');
target_db.is_metabolic = ismember(string(target_db.ENTREZ),recon_model.genes);

for i=1:size(target_db,1)
    entrez= target_db{i,2};
    symbol = dico(find(ismember(dico.ENTREZ,string(entrez))),5);
    target_db(i,3) = unique(table2cell(symbol));
end

%How many metabolic genes for the 24 cancer drugs
target_db_24 = target_db(~ismember(target_db.DrugName,{'Methotrexate','Capecitabine'}),:)
numel(unique(target_db_24(ismember(target_db_24.is_metabolic,1),'ENTREZ'))) %52
numel(unique(string(target_db_24.ENTREZ'))) %52

% Table of druggable pathways
pharma_pathways = table(repmat("NA", numel(genes_ss)*100,1));
pharma_pathways.Gene = repmat("NA", numel(genes_ss)*100,1);
pharma_pathways.Gene_Symbol = repmat("NA", numel(genes_ss)*100,1);

pharma_pathways.total_rxns = [1:numel(genes_ss)*100]';
pharma_pathways.targeted_rxns = [1:numel(genes_ss)*100]';
pharma_pathways.targeted_drugs = [1:numel(genes_ss)*100]';
pharma_pathways.targeted_drugs_names = repmat("NA", numel(genes_ss)*100,1);
pharma_pathways.targeted_rxns_names  = repmat("NA", numel(genes_ss)*100,1);
rxns_associated_genes = recon_model.rxns(~strcmp(recon_model.rules,''));
rxns_associated_genes_idxs = find(~strcmp(recon_model.rules,''));

pathways = string(recon_model.subSystems);
k=1;
for p=1:numel(genes_ss)
    pathway = genes_ss{p};
    %pathway
    pathway_idx = find(ismember(pathways, pathway));
    % select rxns idxs under gene control
    pathway_idx = intersect(rxns_associated_genes_idxs,pathway_idx);
    pathway_rxns = recon_model.rxns(pathway_idx);
    % genes for rxns under gene control
    [~,r_genes_p]=find(recon_model.rxnGeneMat(pathway_idx, :));
    for g=1:numel(r_genes_p)
        gene_idx= r_genes_p(g);
        pharma_pathways.Var1(k) = pathway;
        pharma_pathways.total_rxns(k) = numel(recon_model.rxns(pathway_idx));
        pharma_pathways.Gene(k) =  string(recon_model.genes(gene_idx));
        [gene_rxns_idx ,~]=find(recon_model.rxnGeneMat(:, gene_idx));
        targeted_associated_rxns_idx = intersect(intersect(unique(r_genes),pathway_idx),gene_rxns_idx);
        pharma_pathways.targeted_rxns(k) = numel(targeted_associated_rxns_idx);
        pharma_pathways.targeted_rxns_names(k) = string(strjoin(recon_model.rxns(targeted_associated_rxns_idx),'; '));

        gene_idx_in_target_db = find(ismember(string(target_db.ENTREZ),string(recon_model.genes(gene_idx))));
        if numel(gene_idx_in_target_db)>0
            pharma_pathways.Gene_Symbol(k) =  unique(table2cell(dico(find(ismember(dico.ENTREZ,string(recon_model.genes(gene_idx)))),5)));
            shared_drugs_in_p = unique(target_db(gene_idx_in_target_db,1));
            pharma_pathways.targeted_drugs(k) = numel(shared_drugs_in_p);
            pharma_pathways.targeted_drugs_names(k) = string(strjoin(shared_drugs_in_p.DrugName,'; '));
        end
        k=k+1;
    end
end

writetable(pharma_pathways,'../results/Pharmaceutical_drugs_targets.csv');

% create a table of the number of drugs and rxns for each gene-pathway
%% For herbal drugs
target_db = load('../inputs/Drug_target_all_databases_herbs.mat');
target_db =target_db.drug_target_all_database_herbs;
recon_model = readCbModel('../inputs/Recon3DModel_301.mat');
load('../inputs/dico_201911.mat')
genes = recon_model.genes;
for i=1:numel(genes)
    x = strsplit(genes{i,1},'.');
    genes{i} = x{1};
end
recon_model.genes = genes;
target_genes = intersect(string(unique(target_db.entrez_id)),recon_model.genes);
target_genes_idx =find(ismember(recon_model.genes,target_genes));
recon_model.genes(target_genes_idx)
% all rxns of the all cancer drug targets
[r_genes,~]=find(recon_model.rxnGeneMat(:, target_genes_idx));
numel(unique(r_genes))

% How many unique pathways for the drug targets =45
genes_ss=recon_model.subSystems(unique(r_genes));
genes_ss = unique(string(genes_ss));
numel(unique( genes_ss))

%Mapping gene entrez to symbols
target_db.symbol = string([1:size(target_db.entrez_id,1)]');
target_db.is_metabolic = ismember(string(target_db.entrez_id),recon_model.genes);

for i=1:size(target_db,1)
    entrez= target_db{i,5};
    symbol = dico(find(ismember(dico.ENTREZ,string(entrez))),5);
    target_db(i,6) = unique(table2cell(symbol));
end

% Table of druggable pathways
herbal_pathways = table(repmat("NA", numel(genes_ss)*100,1));
herbal_pathways.Gene = repmat("NA", numel(genes_ss)*100,1);
herbal_pathways.Gene_Symbol = repmat("NA", numel(genes_ss)*100,1);

herbal_pathways.total_rxns = [1:numel(genes_ss)*100]';
herbal_pathways.targeted_rxns = [1:numel(genes_ss)*100]';
herbal_pathways.targeted_drugs = [1:numel(genes_ss)*100]';
herbal_pathways.targeted_drugs_names = repmat("NA", numel(genes_ss)*100,1);
herbal_pathways.targeted_rxns_names  = repmat("NA", numel(genes_ss)*100,1);
rxns_associated_genes = recon_model.rxns(~strcmp(recon_model.rules,''));
rxns_associated_genes_idxs = find(~strcmp(recon_model.rules,''));

pathways = string(recon_model.subSystems);
k=1;
for p=1:numel(genes_ss)
    pathway = genes_ss{p};
    %pathway
    pathway_idx = find(ismember(pathways, pathway));
    % select rxns idxs under gene control
    pathway_idx = intersect(rxns_associated_genes_idxs,pathway_idx);
    pathway_rxns = recon_model.rxns(pathway_idx);
    % genes for rxns under gene control
    [~,r_genes_p]=find(recon_model.rxnGeneMat(pathway_idx, :));
    for g=1:numel(r_genes_p)
        gene_idx= r_genes_p(g);
        herbal_pathways.Var1(k) = pathway;
        herbal_pathways.total_rxns(k) = numel(recon_model.rxns(pathway_idx));
        herbal_pathways.Gene(k) =  string(recon_model.genes(gene_idx));
        [gene_rxns_idx ,~]=find(recon_model.rxnGeneMat(:, gene_idx));
        targeted_associated_rxns_idx = intersect(intersect(unique(r_genes),pathway_idx),gene_rxns_idx);
        herbal_pathways.targeted_rxns(k) = numel(targeted_associated_rxns_idx);
        herbal_pathways.targeted_rxns_names(k) = string(strjoin(recon_model.rxns(targeted_associated_rxns_idx),'; '));

        gene_idx_in_target_db = find(ismember(string(target_db.entrez_id),string(recon_model.genes(gene_idx))));
        if numel(gene_idx_in_target_db)>0
            herbal_pathways.Gene_Symbol(k) =  unique(table2cell(dico(find(ismember(dico.ENTREZ,string(recon_model.genes(gene_idx)))),5)));
            shared_drugs_in_p = unique(target_db(gene_idx_in_target_db,4));
            herbal_pathways.targeted_drugs(k) = numel(shared_drugs_in_p);
            herbal_pathways.targeted_drugs_names(k) = string(strjoin(shared_drugs_in_p.name,'; '));
        end
        k=k+1;
    end
end

writetable(herbal_pathways,'../results/Herbal_drugs_targets.csv');

%% For the 2 tested herbal drugs with target information (Emodin, Scutellarien)
target_db = load('../inputs/Drug_target_all_databases_herbs.mat');
target_db =target_db.drug_target_all_database_herbs;
drugs ={'Emodin'} ;
target_db_tested = target_db(find(ismember(target_db.name,drugs)),:)
recon_model = readCbModel('../inputs/Recon3DModel_301.mat');
load('../inputs/dico_201911.mat')
genes = recon_model.genes;
for i=1:numel(genes)
    x = strsplit(genes{i,1},'.');
    genes{i} = x{1};
end
recon_model.genes = genes;
target_genes = intersect(string(unique(target_db_tested.entrez_id)),recon_model.genes);
target_genes_idx =find(ismember(recon_model.genes,target_genes));
recon_model.genes(target_genes_idx)
% all rxns of the all cancer drug targets
[r_genes,~]=find(recon_model.rxnGeneMat(:, target_genes_idx));
numel(unique(r_genes))

% How many unique pathways for the drug targets =45
genes_ss=recon_model.subSystems(unique(r_genes));
genes_ss = unique(string(genes_ss));
numel(unique( genes_ss))

%Mapping gene entrez to symbols
target_db_tested.symbol = string([1:size(target_db_tested.entrez_id,1)]');
target_db_tested.is_metabolic = ismember(string(target_db_tested.entrez_id),recon_model.genes);

for i=1:size(target_db_tested,1)
    entrez= target_db_tested{i,5};
    symbol = dico(find(ismember(dico.ENTREZ,string(entrez))),5);
    target_db_tested(i,6) = unique(table2cell(symbol));
end

% Table of druggable pathways
herbal_pathways = table(repmat("NA", numel(genes_ss)*100,1));
herbal_pathways.Gene = repmat("NA", numel(genes_ss)*100,1);
herbal_pathways.Gene_Symbol = repmat("NA", numel(genes_ss)*100,1);

herbal_pathways.total_rxns = [1:numel(genes_ss)*100]';
herbal_pathways.targeted_rxns = [1:numel(genes_ss)*100]';
herbal_pathways.targeted_drugs = [1:numel(genes_ss)*100]';
herbal_pathways.targeted_drugs_names = repmat("NA", numel(genes_ss)*100,1);
herbal_pathways.targeted_rxns_names  = repmat("NA", numel(genes_ss)*100,1);
rxns_associated_genes = recon_model.rxns(~strcmp(recon_model.rules,''));
rxns_associated_genes_idxs = find(~strcmp(recon_model.rules,''));

pathways = string(recon_model.subSystems);
k=1;

for p=1:numel(genes_ss)
    pathway = genes_ss{p};
    %pathway
    pathway_idx = find(ismember(pathways, pathway));
    % select rxns idxs under gene control
    pathway_idx = intersect(rxns_associated_genes_idxs,pathway_idx);
    pathway_rxns = recon_model.rxns(pathway_idx);
    % genes for rxns under gene control
    [~,r_genes_p]=find(recon_model.rxnGeneMat(pathway_idx, :));
    for g=1:numel(r_genes_p)
        gene_idx= r_genes_p(g);
        herbal_pathways.Var1(k) = pathway;
        herbal_pathways.total_rxns(k) = numel(recon_model.rxns(pathway_idx));
        herbal_pathways.Gene(k) =  string(recon_model.genes(gene_idx));
        [gene_rxns_idx ,~]=find(recon_model.rxnGeneMat(:, gene_idx));
        targeted_associated_rxns_idx = intersect(intersect(unique(r_genes),pathway_idx),gene_rxns_idx);
        herbal_pathways.targeted_rxns(k) = numel(targeted_associated_rxns_idx);
        herbal_pathways.targeted_rxns_names(k) = string(strjoin(recon_model.rxns(targeted_associated_rxns_idx),'; '));

        gene_idx_in_target_db = find(ismember(string(target_db_tested.entrez_id),string(recon_model.genes(gene_idx))));
        if numel(gene_idx_in_target_db)>0
            herbal_pathways.Gene_Symbol(k) =  unique(table2cell(dico(find(ismember(dico.ENTREZ,string(recon_model.genes(gene_idx)))),5)));
            shared_drugs_in_p = unique(target_db_tested(gene_idx_in_target_db,4));
            herbal_pathways.targeted_drugs(k) = numel(shared_drugs_in_p);
            herbal_pathways.targeted_drugs_names(k) = string(strjoin(shared_drugs_in_p.name,'; '));
        end
        k=k+1;
    end
end

writetable(herbal_pathways,'../results/Emodin_targets.csv');

%% For herbal drugs (Emodin, Scutellarien)
target_db = load('../inputs/Drug_target_all_databases_herbs.mat');
target_db =target_db.drug_target_all_database_herbs;
drugs ={'Scutellarein'} ;
target_db_tested = target_db(find(ismember(target_db.name,drugs)),:)
recon_model = readCbModel('../inputs/Recon3DModel_301.mat');
load('../inputs/dico_201911.mat')
genes = recon_model.genes;
for i=1:numel(genes)
    x = strsplit(genes{i,1},'.');
    genes{i} = x{1};
end
recon_model.genes = genes;
target_genes = intersect(string(unique(target_db_tested.entrez_id)),recon_model.genes);
target_genes_idx =find(ismember(recon_model.genes,target_genes));
recon_model.genes(target_genes_idx)
% all rxns of the all cancer drug targets
[r_genes,~]=find(recon_model.rxnGeneMat(:, target_genes_idx));
numel(unique(r_genes))

% How many unique pathways for the drug targets =45
genes_ss=recon_model.subSystems(unique(r_genes));
genes_ss = unique(string(genes_ss));
numel(unique( genes_ss))

%Mapping gene entrez to symbols
target_db_tested.symbol = string([1:size(target_db_tested.entrez_id,1)]');
target_db_tested.is_metabolic = ismember(string(target_db_tested.entrez_id),recon_model.genes);

for i=1:size(target_db_tested,1)
    entrez= target_db_tested{i,5};
    symbol = dico(find(ismember(dico.ENTREZ,string(entrez))),5);
    target_db_tested(i,6) = unique(table2cell(symbol));
end

% Table of druggable pathways
herbal_pathways = table(repmat("NA", numel(genes_ss)*100,1));
herbal_pathways.Gene = repmat("NA", numel(genes_ss)*100,1);
herbal_pathways.Gene_Symbol = repmat("NA", numel(genes_ss)*100,1);

herbal_pathways.total_rxns = [1:numel(genes_ss)*100]';
herbal_pathways.targeted_rxns = [1:numel(genes_ss)*100]';
herbal_pathways.targeted_drugs = [1:numel(genes_ss)*100]';
herbal_pathways.targeted_drugs_names = repmat("NA", numel(genes_ss)*100,1);
herbal_pathways.targeted_rxns_names  = repmat("NA", numel(genes_ss)*100,1);
rxns_associated_genes = recon_model.rxns(~strcmp(recon_model.rules,''));
rxns_associated_genes_idxs = find(~strcmp(recon_model.rules,''));

pathways = string(recon_model.subSystems);
k=1;

for p=1:numel(genes_ss)
    pathway = genes_ss{p};
    %pathway
    pathway_idx = find(ismember(pathways, pathway));
    % select rxns idxs under gene control
    pathway_idx = intersect(rxns_associated_genes_idxs,pathway_idx);
    pathway_rxns = recon_model.rxns(pathway_idx);
    % genes for rxns under gene control
    [~,r_genes_p]=find(recon_model.rxnGeneMat(pathway_idx, :));
    for g=1:numel(r_genes_p)
        gene_idx= r_genes_p(g);
        herbal_pathways.Var1(k) = pathway;
        herbal_pathways.total_rxns(k) = numel(recon_model.rxns(pathway_idx));
        herbal_pathways.Gene(k) =  string(recon_model.genes(gene_idx));
        [gene_rxns_idx ,~]=find(recon_model.rxnGeneMat(:, gene_idx));
        targeted_associated_rxns_idx = intersect(intersect(unique(r_genes),pathway_idx),gene_rxns_idx);
        herbal_pathways.targeted_rxns(k) = numel(targeted_associated_rxns_idx);
        herbal_pathways.targeted_rxns_names(k) = string(strjoin(recon_model.rxns(targeted_associated_rxns_idx),'; '));

        gene_idx_in_target_db = find(ismember(string(target_db_tested.entrez_id),string(recon_model.genes(gene_idx))));
        if numel(gene_idx_in_target_db)>0
            herbal_pathways.Gene_Symbol(k) =  unique(table2cell(dico(find(ismember(dico.ENTREZ,string(recon_model.genes(gene_idx)))),5)));
            shared_drugs_in_p = unique(target_db_tested(gene_idx_in_target_db,4));
            herbal_pathways.targeted_drugs(k) = numel(shared_drugs_in_p);
            herbal_pathways.targeted_drugs_names(k) = string(strjoin(shared_drugs_in_p.name,'; '));
        end
        k=k+1;
    end
end

writetable(herbal_pathways,'../results/Scutellarein_targets.csv');
