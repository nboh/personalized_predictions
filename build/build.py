import pandas as pd
import micom 


def add_all_probiotic(taxonomy):
    taxonomy_probiotic_total = pd.DataFrame()
    for sample in taxonomy['sample_id'].unique():
        taxonomy_reduced = taxonomy[taxonomy['sample_id'] == sample].copy()
        taxonomy_reduced['abundance'] = taxonomy_reduced['abundance']*0.90
        taxonomy_probiotic = pd.concat([taxonomy_reduced,pd.DataFrame({
                  'sample_id':[sample, sample, sample, sample, sample], 
                  'species': ['Anaerobutyricum hallii', 'Akkermansia muciniphila',
                             'Bifidobacterium longum','Clostridium beijerinckii',
                            'Clostridium butyricum'],
                  'abundance':[0.02, 0.02, 0.02, 0.02, 0.02], 
                  'id':['Anaerobutyricum_hallii', 'Akkermansia_muciniphila',
                             'Bifidobacterium_longum','Clostridium_beijerinckii',
                             'Clostridium_butyricum']})])
        taxonomy_probiotic['sample_id'] = taxonomy_probiotic['sample_id'].astype('str')+'_cocktail'
        taxonomy_probiotic_total = pd.concat([taxonomy_probiotic_total, taxonomy_probiotic])

    return taxonomy_probiotic_total


# locate GEM model database
agora = './agora103_refseq216_species_1.qza'

# convert metagenomics pipeline output to taxonomy table
taxonomy = pd.read_csv('data/S_counts.csv')
taxonomy['sample_id'] = taxonomy['sample_id'].astype('str')
taxonomy.rename(columns = {'fraction_total_reads':'abundance',
                           'name':'species'},
                inplace = True)

taxonomy['id'] = taxonomy['species'].str.replace(' ', '_')

# sum duplicates
taxonomy= taxonomy.groupby(['sample_id','species','id']).sum().reset_index()

# calculate relative abundance
taxonomy['tot_reads'] = taxonomy.groupby('sample_id')['reads'].transform('sum')
taxonomy['abundance'] = taxonomy['reads']/taxonomy['tot_reads']


taxonomy_cocktail = add_all_probiotic(taxonomy)



# combine all models
taxonomy = pd.concat([taxonomy, taxonomy_cocktail])

taxonomy['species'] = taxonomy['species'].str.replace('_',' ')
taxonomy['sample_id'] = taxonomy['sample_id'].astype('str')
taxonomy = taxonomy.groupby(['sample_id','species','id']).sum().reset_index()

# build models with MICOM
manifest = micom.workflows.build(taxonomy, 
                                 model_db=agora, 
                                 out_folder='./model', 
                                 cutoff=0.001, 
                                 threads=10)