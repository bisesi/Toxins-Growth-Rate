import cometspy as c
import cobra
import math
import random
from cobra import Metabolite, Reaction, Model
import pandas as pd
import os
import sys
import itertools
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
import pandas as pd

# define functions
def pick_unique_locations(width, height, n, edge_space = 0):
    locs = []
    while len(locs) < n:
        loc = (random.randrange(edge_space, width - edge_space),
               random.randrange(edge_space, height - edge_space))
        if loc not in locs:
            locs.append(loc)
    return(locs)

def make_RPS_cobra_models_biomass_cost(growth_rate = 1., toxin_cost = 0.02, resistance_cost = 0.0, toxin_prod = 1.):
    # this is different from make RPS cobra models because the toxin and resistance cost is put into the
    # biomass reaction, rather than a separate reaction. this is important because toxins will only
    # be created during growth.also, no need to multiply costs by growth rate. in all cases, they
    # vary directly with growth rate. finally, in all cases as much toxin is produced as carbon is taken
    # up, regardless of growth rate.
    #
    # These toy models use "carbon_c" to grow--that is the only resource
    #
    # the toxin is called "toxin_e"
    '''
    @toxin_prod:  a multiplier.  the multiple of mmol of toxin made per gram of carbon taken up
    '''
    carbon_e = Metabolite(id = "carbon_e",
               compartment = "e")
    carbon_c = Metabolite(id = "carbon_c",
                   compartment = "c")
    EX_carbon_e = Reaction(id = "EX_carbon_e",
                      lower_bound = -growth_rate, # growth rate
                      upper_bound = 1000.)
    EX_carbon_e.add_metabolites({carbon_e: -1})
    carbon_transfer = Reaction(id = "carbon_transfer",
                          lower_bound = 0.,
                          upper_bound = 1000.)
    carbon_transfer.add_metabolites({carbon_e: -1,
                                carbon_c: 1})
    Biomass = Reaction(id = "Biomass",
                  lower_bound = 0.,
                  upper_bound = 1000.)
    Biomass.add_metabolites({carbon_c: -1.})
    # make the toxicity-related metabolites and reactions
    toxin_c = Metabolite(id = "toxin_c", compartment = "c")
    toxin_e = Metabolite(id = "toxin_e", compartment = "e")

    EX_toxin_e = Reaction(id = "EX_toxin_e",
                         lower_bound = -1000.,
                         upper_bound = 1000.)
    EX_toxin_e.add_metabolites({toxin_e: -1})

    toxin_transfer = Reaction(id = "toxin_transfer",
                             lower_bound = -1000.,
                             upper_bound = 0.)
    toxin_transfer.add_metabolites({toxin_e: -1,
                                   toxin_c: 1})
    
    toxin_transfer_sus = Reaction(id = "toxin_transfer_sus",
                             lower_bound = -1000.,
                             upper_bound = 0.)
    toxin_transfer_sus.add_metabolites({toxin_e: 1,
                                   toxin_c: -1})

    Biomass_producer = Reaction(id = "Biomass",
                  lower_bound = 0.,
                  upper_bound = 1000.)
    Biomass_producer.add_metabolites({carbon_c: -(1. + toxin_cost + resistance_cost),
                                     toxin_c: toxin_prod * (1. + toxin_cost + resistance_cost)})

    Biomass_resistant = Reaction(id = "Biomass",
                  lower_bound = 0.,
                  upper_bound = 1000.)
    Biomass_resistant.add_metabolites({carbon_c: -(1. + resistance_cost)})

    producer = Model("producer")
    producer.add_reactions([EX_carbon_e, carbon_transfer, Biomass_producer,
                        EX_toxin_e, toxin_transfer])
    producer.objective = Biomass_producer

    resistant = Model("resistant")
    resistant.add_reactions([EX_carbon_e, carbon_transfer, Biomass_resistant])
    resistant.objective = Biomass_resistant

    susceptible = Model("susceptible")
    susceptible.add_reactions([EX_carbon_e, carbon_transfer, Biomass,
                        EX_toxin_e, toxin_transfer_sus])
    susceptible.objective = Biomass
    return((producer, resistant, susceptible))

# set up treatment parameters
growth_rates = [0.125, 0.25, 0.5, 0.75, 1]
production_costs = [0.01]
toxin_coefficients = [15]
resistance_costs = [0.0]
metabolite_diffs = [5e-6]
toxin_diffs = [5e-7]
add_signal_parameters = ["bounded_linear"]
replicates = [1,2,5]
space_width_starts = [0.01]
cell_density_factor = 2
grid_dim = 50

# put treatments into a dataframe to iterate through
parameter_dict = {'growth_rate': growth_rates, 'production_cost': production_costs, 
                  'toxin_coefficient': toxin_coefficients, 'resistance_cost': resistance_costs,
                  'metabolite_diff': metabolite_diffs, 'toxin_diff': toxin_diffs, 'add_signal_parameter': add_signal_parameters,
                  'replicate': replicates, 'space_widths': space_width_starts}

expanded_grid = pd.DataFrame(itertools.product(*parameter_dict.values()),columns=parameter_dict.keys())

base = f"/Users/abisesi/Desktop/PhD/Projects/Toxins-Growth-Rate/toxin-simulations/simulations_spatial/"

try:
    os.mkdir(base)
except:
    print("Warning: base directory exists. data may be overwritten")

# initialize data frames to store location and final biomass data
total_biomass_data = pd.DataFrame()
metabolite_data = pd.DataFrame()
starting_locations = pd.DataFrame()

# run spatial simulations, storing location and final biomass data before saving as csv after all loops
for i in expanded_grid.index.values:
    growth_rate = expanded_grid.at[i, 'growth_rate']
    production_cost = expanded_grid.at[i, 'production_cost']
    toxin_coefficient = expanded_grid.at[i, 'toxin_coefficient']
    resistance_cost = expanded_grid.at[i, 'resistance_cost']
    metabolite_diff = expanded_grid.at[i, 'metabolite_diff']
    toxin_diff = expanded_grid.at[i, 'toxin_diff']
    add_signal_parameter = expanded_grid.at[i, 'add_signal_parameter']
    replicate = expanded_grid.at[i, 'replicate']
    space_width = expanded_grid.at[i, 'space_widths']
    initial_biomass = 1.e-10 * cell_density_factor
    grid_size = [grid_dim, grid_dim]
    base_hours = 1 

    # iterate through these to figure out which to use, and the spatial seed
    spatial_seed = math.floor(replicate / 1)
    random.seed(spatial_seed) 

    # time step needs to shrink as D increases or dx decreases:
    dt = 0.25 * ((space_width)**2) / metabolite_diff
    # sim hours aka max cycles is hard to pre-determine, you may need to adjust later. 
    hours = base_hours / growth_rate
    max_cycles = int(hours / dt)
    if growth_rate == 0.125 or growth_rate == 0.25:
        max_cycles = 75
    else:
        max_cycles = 30

    producer, resistant, susceptible = make_RPS_cobra_models_biomass_cost(growth_rate = growth_rate, 
                                                                        toxin_cost = production_cost,
                                                                        resistance_cost = resistance_cost,
                                                                        toxin_prod = toxin_coefficient)

    # HERE IS WHERE THE NUMBER OF FOUNDERS IS SET
    locs = pick_unique_locations(grid_size[0], grid_size[1], grid_dim, 3)
    producer_locs = locs[0:5]
    resistant_locs = locs[5:10]
    susceptible_locs = locs[10:]

    producer_models = []
    for l in range(0, len(producer_locs)):
        P = c.model(producer)

        P.obj_style = "MAX_OBJECTIVE_MIN_TOTAL"

        loc = producer_locs[l]
        P.initial_pop = [loc[0], loc[1], initial_biomass]
        producer_models.append(P)

    resistant_models = []
    for k in range(0, len(resistant_locs)):
        R = c.model(resistant)

        R.obj_style = "MAX_OBJECTIVE_MIN_TOTAL"

        ### Fix COMETS thinking biomass reactions with only reactants are exchanges
        R.reactions.loc[R.reactions.REACTION_NAMES == "Biomass", "EXCH"] = False
        index = R.reactions.loc[R.reactions.REACTION_NAMES == "Biomass", "EXCH_IND"].values[0]
        R.reactions.loc[R.reactions.REACTION_NAMES == "Biomass", "EXCH_IND"] = 0
        R.reactions.loc[R.reactions.EXCH_IND > index, "EXCH_IND"] -= 1

        loc = resistant_locs[k]
        R.initial_pop = [loc[0], loc[1], initial_biomass]
        resistant_models.append(R)

    susceptible.add_boundary(susceptible.metabolites.get_by_id("toxin_c"), type="sink")
    susceptible_models = []
    for j in range(0, len(susceptible_locs)):
        S = c.model(susceptible)
        
        S.obj_style = "MAX_OBJECTIVE_MIN_TOTAL"
        
        ### Fix COMETS thinking biomass reactions with only reactants are exchanges
        S.reactions.loc[S.reactions.REACTION_NAMES == "Biomass", "EXCH"] = False
        index = S.reactions.loc[S.reactions.REACTION_NAMES == "Biomass", "EXCH_IND"].values[0]
        S.reactions.loc[S.reactions.REACTION_NAMES == "Biomass", "EXCH_IND"] = 0
        S.reactions.loc[S.reactions.EXCH_IND > index, "EXCH_IND"] -= 1
        
        ### setup the toxicity
        biomass_id = S.reactions.loc[S.reactions.REACTION_NAMES == "Biomass", "ID"].values[0]
        toxin_exch_id = S.reactions.loc[S.reactions.REACTION_NAMES == "EX_toxin_e", "EXCH_IND"].values[0]
        gr_absense_of_toxin = growth_rate
        toxin_conc_where_effect_starts = 0.
        slope_of_toxin_effect = -growth_rate
        toxin_conc_where_effect_saturates = toxin_coefficient
        S.add_signal(biomass_id, toxin_exch_id, 'ub', add_signal_parameter, 
                parms = [gr_absense_of_toxin, 
                        toxin_conc_where_effect_starts, 
                        slope_of_toxin_effect,
                        toxin_conc_where_effect_saturates])

        loc = susceptible_locs[j]
        S.initial_pop = [loc[0], loc[1], initial_biomass]
        susceptible_models.append(S)

    # make params
    p = c.params()
    p.set_param("defaultVmax", 1.)
    p.set_param("defaultKm", 0.000001) 
    p.set_param('maxCycles', max_cycles)
    p.set_param('timeStep', dt)
    p.set_param('spaceWidth', space_width)
    p.set_param('writeMediaLog', True)
    p.set_param('writeBiomassLog', False)
    p.set_param('minSpaceBiomass', 1.e-15)
    p.set_param('defaultDiffConst', metabolite_diff)
    p.set_param('flowDiffRate', 3.e-9)
    p.set_param('writeFluxLog', False)
    p.set_param("totalBiomassLogRate", 1)
    p.set_param("MediaLogRate", 1)
    p.set_param("numDiffPerStep", 1) # this SHOULD be fine since we determined dt based on the faster diff rate

    ## Make the layout
    l = c.layout(susceptible_models + producer_models + resistant_models)
    l.grid = grid_size
    l.set_specific_metabolite("carbon_e", 4.e-7)
    l.media.loc[l.media.metabolite == "toxin_e", "diff_c"] = toxin_diff

    ## run the model, then save total biomass
    sim = c.comets(l, p)
    sim.working_dir = base

    try:
        sim.run(delete_files = False)
    except:
        print(sim.run_output)

    biomass_total = sim.total_biomass
    biomass_total["growth_rate"] = growth_rate
    biomass_total["spatial_seed"] = spatial_seed

    metabolites = sim.media
    metabolites["growth_rate"] = growth_rate
    metabolites["spatial_seed"] = spatial_seed

    total_biomass_data = pd.concat([total_biomass_data, biomass_total])
    metabolite_data = pd.concat([metabolite_data, metabolites])

    spatial_data = pd.DataFrame(columns = ["strain", "x", "y"])
    producer_data = pd.DataFrame({"strain" : "producer",
                                "x" : [x[0] for x in producer_locs],
                                "y" : [x[1] for x in producer_locs]})
    resistant_data = pd.DataFrame({"strain" : "resistant",
                                "x" : [x[0] for x in resistant_locs],
                                "y" : [x[1] for x in resistant_locs]})
    susceptible_data = pd.DataFrame({"strain" : "susceptible",
                                "x" : [x[0] for x in susceptible_locs],
                                "y" : [x[1] for x in susceptible_locs]})
    spatial_data = spatial_data.append(producer_data, 
                                    ignore_index = True).append(resistant_data, 
                                                                ignore_index = True).append(susceptible_data, ignore_index = True)
    spatial_data["growth_rate"] = growth_rate
    spatial_data["spatial_seed"] = spatial_seed

    starting_locations = pd.concat([starting_locations, spatial_data])

    print(i + 1, "of", max(expanded_grid.index.values) + 1)

total_biomass_data.to_csv("/Users/abisesi/Desktop/PhD/Projects/Toxins-Growth-Rate/toxin-simulations/for_reviewers/total_biomass_sink.csv")
metabolite_data.to_csv("/Users/abisesi/Desktop/PhD/Projects/Toxins-Growth-Rate/toxin-simulations/for_reviewers/media_sink.csv")
starting_locations.to_csv("/Users/abisesi/Desktop/PhD/Projects/Toxins-Growth-Rate/toxin-simulations/for_reviewers/starting_locations_sink.csv")

import glob
for f in glob.glob("/Users/abisesi/Desktop/PhD/Projects/Toxins-Growth-Rate/toxin-simulations/for_reviewers/totalbiomasslog*"):
    os.remove(f)

for f in glob.glob("/Users/abisesi/Desktop/PhD/Projects/Toxins-Growth-Rate/toxin-simulations/for_reviewers/*.cmd"):
    os.remove(f)

for f in glob.glob("/Users/abisesi/Desktop/PhD/Projects/Toxins-Growth-Rate/toxin-simulations/for_reviewers/*.txt"):
    os.remove(f)

for f in glob.glob("/Users/abisesi/Desktop/PhD/Projects/Toxins-Growth-Rate/toxin-simulations/for_reviewers/media*"):
    os.remove(f)

for f in glob.glob("/Users/abisesi/Desktop/PhD/Projects/Toxins-Growth-Rate/toxin-simulations/for_reviewers/biomasslog*"):
    os.remove(f)