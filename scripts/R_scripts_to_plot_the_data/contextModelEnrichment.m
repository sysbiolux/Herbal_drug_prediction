function [rxn_per_pathway, P_values] = contextModelEnrichment(contextModel,genericModel)
%contextModelEnrichment A function to calculate the rxn per pathway for
%context specific model

unique_pathways_generic= unique(genericModel.subSystems);
%% calculate the number of reactions in each pathway for generic model
rxs_per_pathway_generic = num2cell(histcounts(categorical(string(genericModel.subSystems))));
count_pathway_generic = [unique_pathways_generic,rxs_per_pathway_generic'];
count_pathway_generic = array2table(count_pathway_generic);
count_pathway_generic.Properties.VariableNames = [{'pathway', 'generic_count'}];

unique_pathways_context = unique(contextModel.subSystems);
rxs_per_pathway_gene_context = num2cell(histcounts(categorical(string(contextModel.subSystems))));
count_pathway_context = [unique_pathways_context,rxs_per_pathway_gene_context'];
count_pathway_context = array2table(count_pathway_context);
count_pathway_context.Properties.VariableNames = [{'pathway', 'context_count'}];
merging_tables = outerjoin(count_pathway_generic, count_pathway_context, 'MergeKeys',true);
%% changing NaN to zero for context_count column
merging_tables.context_count(ismissing(merging_tables.context_count))=0;
%% divide context_count column by generic_count column to obtain rxn_per_pathway percentage
merging_tables.Rxn_per_Pathway =   str2double(merging_tables.context_count) ./ str2double(merging_tables.generic_count);
%% number of total unique reaction in generic model
merging_tables.total_Rxn_generic = repmat(numel(unique(genericModel.rxns)),numel(merging_tables.context_count),1);
%% %% number of total unique reaction in context model
merging_tables.total_Rxn_context = repmat(numel(unique(contextModel.rxns)),numel(merging_tables.context_count),1);

x = merging_tables.context_count;
M = merging_tables.total_Rxn_generic;
K = merging_tables.total_Rxn_context;
N = merging_tables.generic_count;

%merging_tables.P_values = hygecdf(x,M,K,N);
P_values=1;
rxn_per_pathway = merging_tables(:,[1,4]);

end

