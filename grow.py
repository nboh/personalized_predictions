import pandas as pd
import micom

def pectin_supplement(medium):
    medium = pd.concat([medium,
                        pd.DataFrame(
                            index = ['EX_pect_m'],
                            data = {'flux': [1],
                                    'dilution':[1.0],
                                    'metabolite':['pect_m']})])
    return medium

def inulin_supplement(medium):
    medium = pd.concat([medium,
                        pd.DataFrame(
                            index = ['EX_inulin_m'],
                            data = {'flux': [14],
                                    'dilution':[1.0],
                                    'metabolite':['inulin_m']})])
    return medium

def lmn_supplement(medium):
    medium = pd.concat([medium,
                        pd.DataFrame(
                            index = ['EX_inulin_m'],
                            data = {'flux': [14],
                                    'dilution':[1.0],
                                    'metabolite':['lmn30_m']})])
    return medium


def diet_intervention(com_name, diet):
    com = micom.load_pickle(com_name+'.pickle')
    com.medium = diet.flux
    growth = com.cooperative_tradeoff(fraction = 0.7,pfba = False,fluxes = True)
    res = growth.fluxes.mul(growth.members.abundance, axis = 0)
    res = res[res['EX_but(e)']>0]
    sol = res['EX_but(e)'].sum()
    return sol


# load models
com = micom.load_pickle('model/community_baseline.pickle')
com_fp = micom.load_pickle('model/community_fprausnitzii.pickle')
com_am = micom.load_pickle('model/community_amuciniphila.pickle')
com_cb = micom.load_pickle('model/community_cbutyricum.pickle')
com_bi = micom.load_pickle('model/community_binfantis.pickle')

# list of models 
models = [com, com_fp, com_am, com_cb, com_bi]


# load media files
eu_medium = pd.read_csv('diets/eu_medium.csv')
hf_medium = pd.read_csv('diets/hf_medium.csv')

# add fiber interventions
eu_inulin = inulin_supplement(eu_medium)
hf_inulin = inulin_supplement(hf_medium)
eu_pectin = pectin_supplement(eu_medium)
hf_pectin = pectin_supplement(hf_medium)
eu_lmn = lmn_supplement(eu_medium)
hf_lmn = lmn_supplement(hf_medium)

# dietary/prebiotic intervention dict
interventions = {eu_medium:'European',
                 hf_medium:'High-Fiber',
                 eu_inulin:'European + Inulin',
                 hf_inulin:'High-Fiber + Inulin,
                 eu_pectin:'European + Pectin', 
                 hf_pectin:'High-Fiber + Pectin',
                 eu_lmn:'European + Laminarin',
                 hf_lmn:'High-Fiber + Laminarin'}   

# grow combinatorial interventions
res = pd.DataFrame({'community':[],'diet':[],'butyrate':[]})
for community in models:
    for intervention in interventions.keys():
        sol = dietary_intervention(community, interventions[intervention])
        res = pd.concat([res, 
                         pd.DataFrame({'community':[community.name],
                                       'diet':[intervention],
                                       'butyrate':[sol]})
                         
res.to_csv('intervention_results.csv')