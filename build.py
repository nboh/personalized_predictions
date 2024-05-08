import pandas as pd
import micom 

def add_probiotic(taxonomy, probiotic):
    taxonomy['abundance'] = taxonomy['abundance']*0.95
    taxonomy_probiotic = pd.concat([taxonomy,pd.DataFrame({
              'sample_id':[taxonomy['sample_id'].unique()[0]], 
              'genus': [probiotic], 
              'abundance':[0.05], 
              'id':[probiotic.replace(' ','_')]})])
    return taxonomy_probiotic

# locate GEM model database
agora = '/proj/gibbons/refs/micom_databases/v1/agora103_refseq216_species_1.qza'

# convert metagenomics pipeline output to taxonomy table
taxonomy = pd.read_csv('data/S_counts.csv')
taxonomy.rename(columns = {'fraction_total_reads':'abundance',
                           'name':'species', 
                           'sample':'sample_id'},
                inplace = True)
taxonomy['sample_id'] = taxonomy['sample_id'].str.split('-').str[0]
taxonomy['id'] = taxonomy['species'].str.replace(' ', '_')

# sum duplicates
taxonomy= taxonomy.groupby(['sample_id','species','id']).sum().reset_index()

# calculate relative abundance
taxonomy['abundance'] = taxonomy['reads']/taxonomy['reads'].sum(numeric_only=True)

# remove low abundance taxa
taxonomy = taxonomy[taxonomy['abundance'] >=0.001]

# build probiotic taxonomy tables
taxonomy_fp = add_probiotic(taxonomy, 'Faecalibacterium prausnitzii')
taxonomy_cb = add_probiotic(taxonomy, 'Clostridium butyricum')
taxonomy_am = add_probiotic(taxonomy, 'Akkermansia muciniphila')
taxonomy_bi = add_probiotic(taxonomy, 'Bacillus infantis')

# build models with MICOM
com = micom.Community(taxonomy, model_db=agora)
com_fp = micom.Community(taxonomy_fp, model_db=agora)
com_cb = micom.Community(taxonomy_cb, model_db=agora)
com_am = micom.Community(taxonomy_am, model_db=agora)
com_bi = micom.Community(taxonomy_bi, model_db=agora)

# save models
com.to_pickle('model/community_baseline.pickle')
com_fp.to_pickle('model/community_fprausnitzii.pickle')
com_cb.to_pickle('model/community_cbutyricum.pickle')
com_am.to_pickle('model/community_amuciniphila.pickle')
com_bi.to_pickle('model/community_binfantis.pickle')

#