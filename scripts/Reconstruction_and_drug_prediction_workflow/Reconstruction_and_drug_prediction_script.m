clear all
close all
changeCobraSolver('ibm_cplex')

%% TAGS TO CHOOSE WHICH PART OF THE CODE TO RUN
%Run model building tag_Model_Building=1
%Use existing model tag_Model_Building=0
tag_Model_Building =0;
tag_RemoveUnusedGenesFromtheModels=0;
%if 1: repeat the analysis
%if 0:skip the analysis
%if any other value: display the table 1 corresponding of the models
%present in the removedUnusedGenesFolder
%folder
tag_Run_Drug_deletion_herbal_drugs=1;
%if 1: repeat the analysis
%if 0" skip the analysis
%if any other value: display the table with value and the  matlab plot
%To plot the data in R run Plot_Figures_of_the_paper.R

tag_Run_Drug_deletion_cancer_drugs=0;
tag_Run_Drug_deletion_cancer_drugs_on_herbal_models=0;
tag_Run_Drug_deletion_cancer_and_herbal_drugs_on_DMSO_model=0;
tag_Dissimilarity_of_herbal_drugs_to_DMOS=2;
tag_Dissimilarity_of_cancer_drugs_to_DMSO=2;
tag_Similarity_of_herbal_drugs_to_cancer_drugs=2;


%% INPUTS
WORKING_DIR = ''
cd(WORKING_DIR)
cd('./Herbal_drug_prediction/scripts/')

load ('../inputs/metadata.mat');
load('../inputs/colnames.mat');
load('../inputs/Recon3DModel_301.mat');
load ('../inputs/Drug_target_all_databases_herbs.mat');
load '../inputs/GeneTargetPharma.mat'
load '../inputs/Drug_list_cancer_drugs.mat'
load('../inputs/dico_short.mat');
MEM_EBSS_medium = readtable('../inputs/MEM_medium.txt','Delimiter', ' ', 'ReadVariableNames',false);
if exist('../inputs/barcode.txt.gz')>=1 % Extracting the compressed barcode
    gunzip('../inputs/barcode.txt.gz')
end
barcode = readtable('../inputs/barcode.txt', 'ReadRowNames',false); %change ReadRowNames to true if rownames are present
%read row file from barcode
rownames = textread('../inputs/rownames.txt','%s');

List_of_models_to_compare_with={'Capecitabine', 'Methotrexate'};


%% PATHS
%path to Cobratoolbox
%path to FASTCORMICS/rFASTCORMICS
%%
if tag_Model_Building ==1
    
    %extracting condition names (drug names) from metadata
    condition_names = unique(metadata.perturbagen_ch1);
    
    %rename colnames using only the GEO Accession number
    sample_names = colnames;
    for i=1:numel(sample_names)
        name_splitted = strsplit(sample_names{i},'_');
        sample_names(i,1) = cellstr(name_splitted(1,1));
    end
    colnames = sample_names;
    
    
    
    
    %% creation of a new table with ProbeID, EntrezID and index based on rowname file (unmapped genes have '0' as entrez value)
    % create an empty cell array
    table_probe_entrez = table2cell( table(zeros(100000,1),zeros(100000,1),zeros(100000,1),'VariableNames',{'ProbeID', 'ENTREZ','rowname_idx'}));
    
    k = 1;
    
    for t=1:numel(rownames)% for loop for every probe gene ID
        gene_probe=rownames(t); % Determine the probe ID for each t
        if ismember(gene_probe,dico.ProbeID)% Finds if the probe id exists in the dictionary
            gene_entrez_idx = find(ismember(dico.ProbeID,gene_probe));% Determine the row index of the probe id in the dictionary
            gene_entrez = dico.ENTREZ(gene_entrez_idx); % Mapping the correspond entrez id using the probe id index
            %if the Entrez ID is unique fill the table with the probe ID and
            %the corresponding entrez id
            if numel(gene_entrez) ==1
                table_probe_entrez(k,1) = gene_probe;
                table_probe_entrez(k,2) = gene_entrez;
                table_probe_entrez(k,3) = {t};
                k = k +1;
                %if the Entrez ID is not unique, fill the table with the probe ID and
                %the corresponding entrez id repeated
            else
                table_probe_entrez(k:k+numel(gene_entrez)-1,1) = repelem(gene_probe,numel(gene_entrez))';
                table_probe_entrez(k:k+numel(gene_entrez)-1,2) = gene_entrez;
                table_probe_entrez(k:k+numel(gene_entrez)-1,3) = {t};
                k = k + numel(gene_entrez);
                
            end
            %if the Probe ID don't correspond to any Entrez ID , fill the table with the probe ID and 0 for the Entrez id
        else
            table_probe_entrez(k,1) = gene_probe;
            table_probe_entrez(k,2) = {'0'};
            table_probe_entrez(k,3) = {t};
            k = k +1;
        end
    end
    % trasform the cell array to table
    table_probe_entrez= cell2table(table_probe_entrez, 'VariableNames', {'ProbeID', 'Entrez', 'Index'});
    
    % removing the unmapped genes
    % selectig entrez ==0
    indx =  find(ismember(string(table_probe_entrez.Entrez),'0'));
    table_probe_entrez(indx,:)=[];
    
    %Find all duplicated Probe id corresponding to different Entrez ID
    
    % indices to unique values in ProbeID
    [~, ind] = unique(table_probe_entrez.ProbeID);
    % duplicate indices
    duplicate_ind = setdiff(1:size(table_probe_entrez.ProbeID, 1), ind);
    % duplicate values
    duplicate_value = table_probe_entrez.ProbeID(duplicate_ind);
    % find all duplicate probe Id inside table
    all_probe_duplicated =table_probe_entrez.ProbeID(ismember(table_probe_entrez.ProbeID,duplicate_value));
    
    %remove duplicated probe ID and non mapped genes from the table
    table_probe_entrez(ismember(table_probe_entrez.ProbeID,all_probe_duplicated),:) = [];
    
    %Remove duplicated probe ID and non mapped genes from barcode and rownames
    barcode.Properties.RowNames = rownames;
    barcode_filt = barcode(ismember(barcode.Row,table_probe_entrez.ProbeID),:);
    barcode_filt.Properties.RowNames = {};
    barcode.Properties.RowNames = {};
    rownames = table_probe_entrez.ProbeID;
    
    %setting model reconstruction parameters
    %set thresholds
    discretized = table2array(barcode_filt);
    
    expression_threshold      = 5; %5 frma
    inexpression_threshold    = 0; %0 frma
    
    match_expressed     = table2array(barcode_filt) > expression_threshold;
    match_unexpressed   = table2array(barcode_filt) <= inexpression_threshold;
    match_unknown       = table2array(barcode_filt) <= expression_threshold & table2array(barcode_filt) > inexpression_threshold;
    %apply thresholds
    discretized(match_expressed)   = 1;
    discretized(match_unexpressed) = -1;
    discretized(match_unknown)     = 0;
    
    figure
    imagesc(discretized)
    clear data expression_threshold inexpression_threshold match_expressed match_unexpressed match_unknown
    
    
    
    
    %Setup model
    Recon3_model = Recon3DModel; clear Recon3DModel
    Recon3_model.rev=zeros(numel(Recon3_model.lb),1);
    Recon3_model.rev(Recon3_model.lb<0)=1;
    %Remove gene version from the model
    Recon3_model.genes=regexprep(Recon3_model.genes,'\.[0-9]+$','');
    
    epsilon = 1e-4;
    consensus_proportion = 0.9; %gene has to be expressed in 90% of the cases in order to be included.
    % Only relevant if you want to create one model from different samples
    already_mapped_tag = 0;
    
    %read the table with medium components
    MEM_EBSS_medium = table2array(MEM_EBSS_medium);
    
    unpenalizedSystems = {'Transport, endoplasmic reticular';
        'Transport, extracellular';
        'Transport, golgi apparatus';
        'Transport, mitochondrial';
        'Transport, peroxisomal';
        'Transport, lysosomal';
        'Transport, nuclear'};
    %unique(string(Recon3_model.subSystems))
    unpenalized = Recon3_model.rxns(ismember((vertcat(Recon3_model.subSystems {:})),unpenalizedSystems));
    
    optional_settings.unpenalized = unpenalized;
    optional_settings.func = {'biomass_reaction'};
    biomass_rxn = 'biomass_reaction';
    
    optional_settings.medium = MEM_EBSS_medium;
    
    %mkdir models_
    % Create generic models
    % create models_keep_generic
    %dico = cellstr(table2array(dico));
    dico = cellstr(table2array(table_probe_entrez(:,1:2)));
    
    %changing the Four mixture name in simply 'Mixture'
    strg_to_replace = {'Four mixture(Tanshinone IIA:Salvianic acid A sodium:Protocatechuic aldehyde:Salvianolic acid B=1:1:1:1)', '(+)2-(1-hydroxyl-4-oxocyclohexyl) ethyl caffeate'};
    new_names = {'Four mixture', 'Ethyl caffeate'};
    condition_names = replace(condition_names, strg_to_replace, new_names) ;
    metadata.perturbagen_ch1 = replace(metadata.perturbagen_ch1, strg_to_replace, new_names) ;
    
    models_keep_generic = zeros(numel(Recon3_model.rxns), numel(condition_names));
    for i=1:numel(condition_names)
        clear context_model
        %define the conditions
        drug = condition_names(i);
        
        % define the GSE_Id to the corresponding drug inside metadata
        GSE_for_i_conditions = metadata.geo_accession(ismember(metadata.perturbagen_ch1, drug));
        if isempty(GSE_for_i_conditions)
            sss
        end
        % define the GSE index of the barcode inside colnames
        GSE_idx_barcode = find(ismember(colnames,GSE_for_i_conditions));
        [context_model, A_final] = fastcormics(Recon3_model, discretized(:,GSE_idx_barcode), rownames, dico, biomass_rxn, already_mapped_tag, consensus_proportion, epsilon, optional_settings);
        context_model = removeUnusedGenes(context_model);
        name=strcat('../models/original models/',condition_names(i),'_model.mat');
        save(name{1},'context_model');
        
    end
    
    % for i=1:numel(condition_names)
    %     clear context_model
    %     load('./models/'+string(condition_names{i})+'_model.mat','context_model');
    %     A_final = find(ismember(Recon3_model.rxns,context_model.rxns));
    %     models_keep_generic(A_final,i) = 1;
    % end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part 2 Numerics of the models (Table 1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
condition_names = unique(metadata.perturbagen_ch1);
strg_to_replace = {'Four mixture(Tanshinone IIA:Salvianic acid A sodium:Protocatechuic aldehyde:Salvianolic acid B=1:1:1:1)', '(+)2-(1-hydroxyl-4-oxocyclohexyl) ethyl caffeate'};
new_names = {'Four mixture', 'Ethyl caffeate'};
condition_names = replace(condition_names, strg_to_replace, new_names) ;
Prediction_summary=zeros(numel(condition_names),4);
condition_names_T = condition_names';
save('../results/condition_names.mat','condition_names_T', '-V6');

if tag_RemoveUnusedGenesFromtheModels==1
    
    model_tested = struct();
    A_keep=zeros(numel(Recon3DModel.rxns), numel(condition_names));
    Number_of_mets=zeros(numel(condition_names),1);
    Number_of_rxns=zeros(numel(condition_names),1);
    Number_of_genes=zeros(numel(condition_names),1);
    
    for i=1:numel(condition_names)
        name=strcat('../models/original models/',condition_names(i),'_model.mat');
        model=load(name{1});
        model=removeUnusedGenes(model.context_model);
        model_tested(i).model=model;
        A_keep(ismember(Recon3DModel.rxns,model.rxns),i)=1;
        Number_of_mets(i)= numel(model.mets);
        Number_of_rxns(i)= numel(model.rxns);
        Number_of_genes(i)=numel(model.genes);
        name=strcat('../models/models after removal of unused genes/',condition_names(i),'_removed_unused_genes_model.mat');
        save(name{1},'model')
        
    end
    
    clear ans new_names i strg_to_replace name model_tested model
    save ('../results/Presence_Matrix','A_keep')
elseif tag_RemoveUnusedGenesFromtheModels==0
    load('../results/Presence_Matrix','A_keep')
elseif tag_RemoveUnusedGenesFromtheModels~=0 &&  tag_RemoveUnusedGenesFromtheModels~=1
    model_tested = struct();
    A_keep=zeros(numel(Recon3DModel.rxns), numel(condition_names));
    Number_of_mets=zeros(numel(condition_names),1);
    Number_of_rxns=zeros(numel(condition_names),1);
    Number_of_genes=zeros(numel(condition_names),1);
    for i=1:numel(condition_names)
        name=strcat('../models/models after removal of unused genes/',condition_names(i),'_removed_unused_genes_model.mat');
        load(name{1});
        model_tested(i).model=model;
        clear model
        
        A_keep(ismember(Recon3DModel.rxns,model.rxns),i)=1;
        Number_of_mets(i)= numel(model.mets);
        Number_of_rxns(i)= numel(model.rxns);
        Number_of_genes(i)=numel(model.genes);
        save('../results/Presence_Matrix','A_keep')
    end
    
    Numerics(1,1)=mean(Number_of_rxns);
    Numerics(2,1)=min(Number_of_rxns);
    Numerics(3,1)=max(Number_of_rxns);
    
    Numerics(1,2)=mean(Number_of_mets);
    Numerics(2,2)=min(Number_of_mets);
    Numerics(3,2)=max(Number_of_mets);
    
    Numerics(1,3)=mean(Number_of_genes);
    Numerics(2,3)=min(Number_of_genes);
    Numerics(3,3)=max(Number_of_genes);
    Numerics=table(Numerics(:,1), Numerics(:,2),Numerics(:,3));
    Numerics.Properties.VariableNames={'reactions', 'metabolites', 'genes'};
    Numerics;
end
clear strg_to_replace i name model new_names model_tested


if tag_Run_Drug_deletion_herbal_drugs==1
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Part 3 Drug Deletion (Table 1)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    match=contains(condition_names,'DMSO');
    if ~exist('DMSO_model','var')
        name=strcat('../models/models after removal of unused genes/',condition_names(match),'_removed_unused_genes_model.mat');
        load(name{1});
        DMSO_model=model;
        clear model
    end
    DrugList_herbal=unique(drug_target_all_database_herbs.name);
    DMSO_model.c=zeros(numel(DMSO_model.rxns),1);
    DMSO_model.c(ismember(DMSO_model.rxns, 'biomass_reaction'))=1;
    grR_keep=ones(numel(DrugList_herbal),1)*-1;
    solwt=optimizeCbModel(DMSO_model);
    for i=1: numel(DrugList_herbal)
        geneID=drug_target_all_database_herbs.entrez_id(ismember(drug_target_all_database_herbs.name, DrugList_herbal(i)));
        geneID_cell=cell(size(geneID));
        for j=1:numel(geneID)
            geneID_cell{j}=num2str(geneID(j));
        end
        geneID_cell= intersect(DMSO_model.genes,geneID_cell);
        [modelDel,~,constrRxnNames] = deleteModelGenes(DMSO_model,geneID_cell);
        
        solmt=optimizeCbModel(modelDel);
        grR=solmt.f/solwt.f;
        grR_keep(i)=grR;
    end
    
    KO_herbal=DrugList_herbal(grR_keep<0.99);
    clear i grR_keep grR modelDel
    Prediction_summary(ismember(condition_names, KO_herbal),1)=1;
    
end


%% KO breast cancer on DMSO
if tag_Run_Drug_deletion_cancer_drugs==1
    if ~exist('DMSO_model','var')
        match=contains(condition_names,'DMSO');
        name=strcat('../models/models after removal of unused genes/',condition_names(match),'_removed_unused_genes_model.mat');
        load(name{1});
        DMSO_model=model;
        clear model
    end
    DMSO_model.c=zeros(numel(DMSO_model.rxns),1);
    DMSO_model.c(ismember(DMSO_model.rxns, 'biomass_reaction'))=1;
    grR_keep=ones(numel(Drug_list_pharmaceutical),1)*-1;
    solwt=optimizeCbModel(DMSO_model);
    
    
    for i=1: numel(Drug_list_pharmaceutical)
        geneID=cellstr(GeneTargetPharma.ENTREZ(ismember(GeneTargetPharma.DrugName, Drug_list_pharmaceutical(i))));
        geneID=intersect(DMSO_model.genes, geneID);
        [modelDel,~,constrRxnNames] = deleteModelGenes(DMSO_model,geneID);
        
        solmt=optimizeCbModel(modelDel);
        grR=solmt.f/solwt.f;
        grR_keep(i)=grR;
    end
    
    KO_pharma=Drug_list_pharmaceutical(grR_keep<0.99);
    
    clear i grR_keep grR
end
%% Cancer Drugs on herbal models
if tag_Run_Drug_deletion_cancer_drugs_on_herbal_models==1
    
    grRatio_keep=ones(numel(Drug_list_pharmaceutical),numel(condition_names))*1;
    grRateWT_keep=ones(numel(Drug_list_pharmaceutical),1);
    
    for i=1:numel(condition_names)
        name=strcat('../models/models after removal of unused genes/',condition_names(i),'_removed_unused_genes_model.mat');
        load(name{1});
        model_tmp=model;
        clear model
        model_tmp.c=zeros(numel(model_tmp.rxns),1);
        model_tmp.c(ismember(model_tmp.rxns, 'biomass_reaction'))=1;
        %setting biomass reaction as objective function
        % model optimizagion
        % call DrugDeletion function
        [grRatio, grRateKO, grRateWT, ~, ~, ~] = DrugDeletion_v2(model_tmp, 'FBA', Drug_list_pharmaceutical, GeneTargetPharma );
        grRateWT_keep(i)=grRateWT;
        grRatio_keep(:,i)=grRatio;
        
    end
    
    
    [r,c]=find(grRatio_keep<0.99);
    combinations=[r,c];
    clear i grRatio_keep grRatio c r
    
end
% No other than Capecitabine and Meso

%%
if tag_Run_Drug_deletion_cancer_and_herbal_drugs_on_DMSO_model==1
    
    grR_keep=ones(size(DMSO_model.rxns,1),numel(condition_names))*1;
    
    for i=1:numel(Drug_list_pharmaceutical)
        for j=1: numel(DrugList_herbal)
            geneID=cellstr(GeneTargetPharma.ENTREZ(ismember(GeneTargetPharma.DrugName, Drug_list_pharmaceutical(i))));
            geneID=intersect(DMSO_model.genes, geneID);
            geneID2=drug_target_all_database_herbs.entrez_id(ismember(drug_target_all_database_herbs.name, DrugList_herbal(j)));
            geneID=union(geneID, cellstr(string(geneID2)));
            geneID=intersect(DMSO_model.genes, geneID);
            [modelDel,~,constrRxnNames] = deleteModelGenes(DMSO_model,geneID);
            
            solmt=optimizeCbModel(modelDel);
            grR=solmt.f/solwt.f;
            grR_keep(i,j)=grR;
        end
        
    end
    
    
    % remove combination known from the KO_drugs
    grR_keep(ismember(Drug_list_pharmaceutical, KO_pharma),:)=1;
    % remove combination known from the KO_drugs
    grR_keep(:,ismember(DrugList_herbal, KO_herbal))=1;
    [r,c]=find(grR_keep<0.99);
    
    combinations2=[r,c];
    clear c r i j
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part 4 DISSIMILARITY OF THE HERBAL TO DMSO (Figure 1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if tag_Dissimilarity_of_herbal_drugs_to_DMOS~=0
    
    %a)herbal models
    
    J2 = squareform(pdist(A_keep','jaccard'));
    Dissimilarity_score=table(condition_names,J2(:,contains(condition_names, 'DMSO')));
    
    Dissimilarity_score = sortrows(Dissimilarity_score,'Var2','ascend');
 figure(1)
    barh(table2array(Dissimilarity_score(end-10:end,2)))
    yticklabels(Dissimilarity_score.condition_names(end-10:end,1))

    Dissimilarity_score = sortrows(Dissimilarity_score,'Var2','descend');

    
    clear J
    Prediction_summary(ismember(condition_names, table2array(Dissimilarity_score(1:5,1))),2)=1;
    save('../results/Herbal_Dissimilarity_to_DMSO.mat','Dissimilarity_score')
    
end







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part 4 DISSIMILARITY OF THE CANCER MODELS TO DMSO (Figure 2a)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if tag_Dissimilarity_of_cancer_drugs_to_DMSO==1
    
    %cancer drugs models
    disp('FVA cancer drugs')
    FVA_results_cancer_drugs=struct();
    for i=1:numel(Drug_list_pharmaceutical)
        disp(strcat(num2str(i),'/', num2str(numel(Drug_list_pharmaceutical))));
        model= DMSO_model;
        geneID=cellstr(GeneTargetPharma.ENTREZ(ismember(GeneTargetPharma.DrugName, Drug_list_pharmaceutical(i))));
        geneID=intersect(DMSO_model.genes, geneID);
        [modelDel,~,constrRxnNames] = deleteModelGenes(model,geneID);
        
        [minFluxtest,maxFluxtest] = fluxVariability(modelDel,100);
        FVA_results_cancer_drugs(i).min=minFluxtest;
        FVA_results_cancer_drugs(i).max=maxFluxtest;
        FVA_results_cancer_drugs(i).rxns= model.rxns;
    end
    save('../results/FVA_cancer_drugs.mat','FVA_results_cancer_drugs')
end

disp('FVA herbal drugs')
if tag_Dissimilarity_of_herbal_drugs_to_DMOS==1
    
    
    FVA_results_herbal_drugs = struct();
    
    for i=1:numel(condition_names)
        disp(strcat(num2str(i),'/', num2str(numel(condition_names))));
        name=strcat('../models/models after removal of unused genes/',condition_names(i),'_removed_unused_genes_model.mat');
        load(name{1});
        model_tested=model;
        clear model
        model_tested.c=zeros(numel(model_tested.rxns),1);
        model_tested.c(ismember(model_tested.rxns, 'biomass_reaction'))=1;
        [minFluxtest,maxFluxtest] = fluxVariability(model_tested,100);
        FVA_results_herbal_drugs(i).min=minFluxtest;
        FVA_results_herbal_drugs(i).max=maxFluxtest;
        FVA_results_herbal_drugs(i).rxns=model_tested.rxns;
        
    end
    save('../results/FVA_herbal_drugs.mat','FVA_results_herbal_drugs')
else
    load('../results/FVA_herbal_drugs.mat','FVA_results_herbal_drugs')

end
if tag_Dissimilarity_of_cancer_drugs_to_DMSO~=0
    res_keep=zeros(numel(Drug_list_pharmaceutical),1);
    
    if tag_Dissimilarity_of_cancer_drugs_to_DMSO~=1
        load('../results/FVA_cancer_drugs.mat','FVA_results_cancer_drugs')
        load('../results/FVA_herbal_drugs.mat','FVA_results_herbal_drugs')
        load('../models/models after removal of unused genes/DMSO_removed_unused_genes_model.mat')
        DMSO_model=model;
    end
    for i=1:numel(Drug_list_pharmaceutical)
        all_reactions=unique([FVA_results_cancer_drugs(i).rxns;DMSO_model.rxns]);
        v= zeros(numel(all_reactions),4);
        v1mins=v(:,1);
        v1maxs=v(:,2);
        v2mins=v(:,3);
        v2maxs=v(:,4);
        [~, IA,IB]=intersect(all_reactions,FVA_results_cancer_drugs(i).rxns);
        
        v1mins(IA)=FVA_results_cancer_drugs(i).min(IB);
        v1maxs(IA)=FVA_results_cancer_drugs(i).max(IB);
        match=contains(condition_names, 'DMSO');
        
        [I, IA,IB]=intersect(all_reactions,FVA_results_herbal_drugs(match).rxns);
        
        v2mins(IA)=FVA_results_herbal_drugs(match).min(IB);
        v2maxs(IA)=FVA_results_herbal_drugs(match).max(IB);
        
        res = FVA_similarity_Thomas(v1mins, v1maxs, v2mins, v2maxs);
        res_keep(i)=res;
        
        
    end
    
    figure(2)
    
    Dissimilarity_to_DMSO=table(Drug_list_pharmaceutical,1-res_keep);
    Dissimilarity_to_DMSO = sortrows(Dissimilarity_to_DMSO,'Var2','ascend');
    barh(table2array(Dissimilarity_to_DMSO(:,2)))
    xlabel('Dissimilarity between the DMSO models constrained with breast cancer drugs and the unconstrained models')
    ylabel('Breast cancer drugs')
    save('../results/Pharma_Dissimilarity_to_DMSO.mat','Dissimilarity_to_DMSO')
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Part 5 SIMILARITY TO CANCER DRUGS (Figure 2b)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if  tag_Similarity_of_herbal_drugs_to_cancer_drugs~=0
    if tag_Similarity_of_herbal_drugs_to_cancer_drugs~=1
        load('../results/FVA_cancer_drugs.mat','FVA_results_cancer_drugs')
        load('../results/FVA_herbal_drugs.mat','FVA_results_herbal_drugs')
    end
    Similarity=struct();
    t_keep=[];
    for j=1:numel(List_of_models_to_compare_with)
        
        index_model_2_compare= contains( Drug_list_pharmaceutical,List_of_models_to_compare_with{j});
        res_keep=zeros(numel(condition_names),1);
        
        for i=1:numel(condition_names)
            all_reactions=unique([FVA_results_herbal_drugs(i).rxns;FVA_results_cancer_drugs(index_model_2_compare).rxns]);
            v= zeros(numel(all_reactions),4);
            v1mins=v(:,1);
            v1maxs=v(:,2);
            v2mins=v(:,3);
            v2maxs=v(:,4);
            [~, IA,IB]=intersect(all_reactions,FVA_results_herbal_drugs(i).rxns);
            
            v1mins(IA)=FVA_results_herbal_drugs(i).min(IB);
            v1maxs(IA)=FVA_results_herbal_drugs(i).max(IB);
            [~, IA,IB]=intersect(all_reactions,FVA_results_cancer_drugs(index_model_2_compare).rxns);
            v2mins(IA)=FVA_results_cancer_drugs(index_model_2_compare).min(IB);
            v2maxs(IA)=FVA_results_cancer_drugs(index_model_2_compare).max(IB);
            
            res = FVA_similarity_Thomas(v1mins, v1maxs, v2mins, v2maxs);
            res_keep(i)=res;
            
            
        end
        t=table(condition_names, res_keep);
        Similarity(j).T=t;
        t = sortrows(t,'res_keep','descend');
        t_keep=[t_keep;t(1:10,1)];
        t_keep=unique(t_keep);
        
        Prediction_summary(ismember(condition_names,table2array(t_keep)),3)=1;
    end
    Similarity_to_cancer_drugs= Similarity(1).T;
    
    for i=2:numel(Similarity)
        [~,IA, IB]= intersect(Similarity_to_cancer_drugs.condition_names, Similarity(i).T.condition_names);
        Similarity_to_cancer_drugs(IA,i+1)=table(Similarity(i).T.res_keep(IB));
        
    end
    for i=1:numel(Similarity)
        name=strcat('../results/Similarity_',List_of_models_to_compare_with(i));
        writetable(Similarity(i).T, name{1})
    end
    
    Similarity_to_cancer_drugs(:,end+1)=table(mean(table2array(Similarity_to_cancer_drugs(:,2:end)),2));
    Similarity_to_cancer_drugs = sortrows(Similarity_to_cancer_drugs,'Var4','descend');
    
    
   
    match=ismember(Similarity_to_cancer_drugs(:,1),t_keep);
    cgo=clustergram(table2array(Similarity_to_cancer_drugs(match,2:3)), 'RowLabels',table2array(Similarity_to_cancer_drugs(match,1)),...
    'ColumnLabels',List_of_models_to_compare_with ,...
    'ColumnLabelsRotate',270, ...
    'DisplayRange',30,...
    'Cluster', 'ROW', ...
    'symmetric','False', 'Colormap','winter');  
    plot(cgo)

    save('../results/Similarity_to_cancer_drugs.mat','Similarity_to_cancer_drugs')
    
    clear res_keep res t t_keep  i IA IB I v1maxs v2
    
end



T_Prediction_summary=table(condition_names,Prediction_summary);

T_Prediction_summary=T_Prediction_summary(sum(Prediction_summary,2)~=0,:);

save('../results/Prediction_summary','T_Prediction_summary')



