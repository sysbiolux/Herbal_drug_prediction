%% Determine the reaction presence for 3 tested herbal drugs + dmso in androgen pathway
WORKING_DIR = ''
cd(WORKING_DIR)
cd('./Herbal_drug_prediction/scripts/')

brc_model = readCbModel('../models/original models/Bruceine D_model.mat');
emd_model = readCbModel('../models/original models/Emodin_model.mat');
sct_model  = readCbModel('../models/original models/Scutellarein_model.mat');
dmso_model  = readCbModel('../models/original models/DMSO_model.mat');
recon_model = readCbModel('../inputs/Recon3DModel_301.mat');

columns = {'Pathway','RXN_id','Bruceine_D','Emodin','Scutellarein','DMSO'};
% pathways that are targetted by most cancer drugs
cancer_pathways = readtable('../results/Summary_Pharmaceutical_drugs_targets.csv');
cancer_pathways_list = cancer_pathways{cancer_pathways.targeted_drugs>15,'genes_ss'};
herbal_rxns_Table = table();
herbal_rxns_Table.Properties.VariableNames;
pathway={'Androgen and estrogen synthesis and metabolism'};
for p=1:numel(cancer_pathways_list)
    pathway = cancer_pathways_list(p);

    brc_rxns = brc_model.rxns(find(ismember(string(brc_model.subSystems),pathway)));
    emd_rxns = emd_model.rxns(find(ismember(string(emd_model.subSystems),pathway)));
    sct_rxns = sct_model.rxns(find(ismember(string(sct_model.subSystems),pathway)));
    dmso_rxns = dmso_model.rxns(find(ismember(string(dmso_model.subSystems),pathway)));
    recon_rxns = recon_model.rxns(find(ismember(string(recon_model.subSystems),pathway)));
    RXN_id = unique([brc_rxns' emd_rxns' sct_rxns' dmso_rxns'])';

    herbal_rxns_T = table(RXN_id);%repmat("NA", numel(genes_ss),1));
    herbal_rxns_T.Pathway = repmat(pathway,numel(RXN_id),1);
    herbal_rxns_T.Bruceine_D = ismember(RXN_id,brc_rxns);
    herbal_rxns_T.Emodin = ismember(RXN_id,emd_rxns);
    herbal_rxns_T.Scutellarein = ismember(RXN_id,sct_rxns);
    herbal_rxns_T.DMSO = ismember(RXN_id,dmso_rxns);
    herbal_rxns_Table = [herbal_rxns_Table;herbal_rxns_T];
end
writetable(herbal_rxns_Table,'../results/Tested_drugs_rxn_presence.csv');