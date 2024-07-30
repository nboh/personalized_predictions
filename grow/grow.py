import pandas as pd
import micom
import micom.measures
import itertools

               
               
fiber_flux  = {
'maltodextrin':{'dextrin': 70},
'pectin':{'pect': 1},
'hemp':{'cellul': 0.845},
'HMOs':{'fuc_L':140, 
       'gal':140,
       'glc_D':140},
'polydextrose':{'glc_D': 422.5},
'inulin':{'inulin': 14},
'starch':{'strch1': 38.5},
'acacia':{'arab_L':150, 
       'gal':150,
       'rmn':150}
}


def add_fiber_supplements(medium, fiber_flux):
    media = {}
    media['baseline'] = medium.copy()
    
    # Generate all combinations of 1 to 8 fibers
    fibers = list(fiber_flux.keys())
    for r in range(1, len(fibers) + 1):
        for combination in itertools.combinations(fibers, r):
            # Create a label for this combination
            combination_label = '+'.join(combination)
            
            # Start with a copy of the original medium
            medium_temp = medium.copy()
            
            # Add each fiber in the combination
            for fiber in combination:
                compounds = {k: v / r for k, v in fiber_flux[fiber].items()}
                compounds_df = pd.DataFrame(compounds.items(), columns=['metabolite', 'flux'])
                compounds_df['reaction'] = 'EX_' + compounds_df['metabolite'] + '_m'
                compounds_df['metabolite'] = compounds_df['metabolite'] + '_m'
                medium_temp = pd.concat([medium_temp, compounds_df])
            
            # Group by reaction and metabolite, summing fluxes
            medium_temp = medium_temp.groupby(['reaction', 'metabolite']).sum(numeric_only=True).reset_index()
            media[combination_label] = medium_temp
    
    return media

# load model manifest
manifest = pd.read_csv('model_pendulum/manifest.csv')


# load media files
eu_medium = pd.read_csv('diets/european_medium.csv')
med_medium = pd.read_csv('diets/mediterranean_diet_v2.csv')
bad_medium = pd.read_csv('diets/unhealthy_diet_v2.csv')


# add fiber interventions
eu_media = add_fiber_supplements(eu_medium, fiber_flux)
med_media = add_fiber_supplements(med_medium, fiber_flux)
bad_media = add_fiber_supplements(bad_medium, fiber_flux)


# dietary/prebiotic intervention dict
diets = {
         'European':eu_media,
         'Mediterranean':med_media,
         'Unhealthy':bad_media
        }

probiotics = ['Akkermansia_muciniphila','Anaerobutyricum_hallii', 'Clostridium_butyricum','Clostridium_beijerinckii', 'Bifidobacterium_longum']

# grow combinatorial interventions
res_but = pd.DataFrame()
res_rates = pd.DataFrame()
for diet in diets.keys():
    media = diets[diet]
    for intervention in media.keys():
        growth = micom.workflows.grow(manifest, 
                                      model_folder = './model_pendulum/',
                                      medium = media[intervention],
                                      tradeoff = 0.95, 
                                      strategy = 'none',
                                      threads = 20)
        production = micom.measures.production_rates(growth)  
        
        but = production[production['name'] == 'butyrate']
        but.loc[:,'prebiotic'] = intervention
        but.loc[:,'diet'] = diet
        rates = growth.growth_rates[growth.growth_rates.index.isin(probiotics)]
        rates.loc[:,'prebiotic'] = intervention
        rates.loc[:,'diet'] = diet
        res_but = pd.concat([res_but, but])
        res_rates = pd.concat([res_rates, rates])
        res_but.to_csv('butyratePredictions_v2.csv')
        res_rates.to_csv('growthRates_v2.csv')