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
                        EX_toxin_e])
    susceptible.objective = Biomass
    return((producer, resistant, susceptible))

## set up treatment parameters
growth_rates = [0.125, 0.25, 0.5, 0.75, 1]
production_costs = [0.01]
toxin_coefficients = [15]
resistance_costs = [0.0]
metabolite_diffs = [5e-6]
toxin_diffs = [5e-7]
add_signal_parameters = ["bounded_linear"] #linear, generalized_logistic
space_width_starts = [0.01]
cell_density_factor = 2
grid_dim = 50

# put treatments into a dataframe to iterate through
parameter_dict = {'growth_rate': growth_rates, 'production_cost': production_costs, 
                  'toxin_coefficient': toxin_coefficients, 'resistance_cost': resistance_costs,
                  'metabolite_diff': metabolite_diffs, 'toxin_diff': toxin_diffs, 'add_signal_parameter': add_signal_parameters,
                  'space_widths': space_width_starts}

expanded_grid = pd.DataFrame(itertools.product(*parameter_dict.values()),columns=parameter_dict.keys())

base = f"/Users/abisesi/Desktop/PhD/Projects/Toxins-Growth-Rate/toxin-simulations/simulations_liquid/"

try:
    os.mkdir(base)
except:
    print("Warning: base directory exists. data may be overwritten")

# initialize data frames to store final biomass data
liquid_sim_data_final = pd.DataFrame()

# run liquid simulations, storing final biomass data before saving as csv after all loops
for i in expanded_grid.index.values:
    growth_rate = expanded_grid.at[i, 'growth_rate']
    production_cost = expanded_grid.at[i, 'production_cost']
    toxin_coefficient = expanded_grid.at[i, 'toxin_coefficient']
    resistance_cost = expanded_grid.at[i, 'resistance_cost']
    metabolite_diff = expanded_grid.at[i, 'metabolite_diff']
    add_signal_parameter = expanded_grid.at[i, 'add_signal_parameter']
    space_width = expanded_grid.at[i, 'space_widths']
    initial_biomass = 1.e-10 * cell_density_factor
    toxin_diff = expanded_grid.at[i, 'toxin_diff']
    base_hours = 1 

    if growth_rate in [0.125, 0.25]:
        maxCycles = 500
    else:
        maxCycles = 300

    producer, resistant, susceptible = make_RPS_cobra_models_biomass_cost(growth_rate = growth_rate, 
                                                                        toxin_cost = production_cost,
                                                                        resistance_cost = resistance_cost,
                                                                        toxin_prod = toxin_coefficient)

    P = c.model(producer)
    R = c.model(resistant)
    S = c.model(susceptible)

    P.obj_style = "MAX_OBJECTIVE_MIN_TOTAL"
    R.obj_style = "MAX_OBJECTIVE_MIN_TOTAL"
    S.obj_style = "MAX_OBJECTIVE_MIN_TOTAL"

    ### Fix COMETS thinking biomass reactions with only reactants are exchanges
    R.reactions.loc[R.reactions.REACTION_NAMES == "Biomass", "EXCH"] = False
    index = R.reactions.loc[R.reactions.REACTION_NAMES == "Biomass", "EXCH_IND"].values[0]
    R.reactions.loc[R.reactions.REACTION_NAMES == "Biomass", "EXCH_IND"] = 0
    R.reactions.loc[R.reactions.EXCH_IND > index, "EXCH_IND"] -= 1

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

    S.initial_pop = [0, 0, 8e-9]
    R.initial_pop = [0, 0, 1e-9]
    P.initial_pop = [0, 0, 1e-9]

    l = c.layout([P, R, S])
    l.set_specific_metabolite("carbon_e", 4.e-7)

    p = c.params()
    p.set_param("defaultVmax", 1.)
    p.set_param("defaultKm", 0.000001) 
    p.set_param('maxCycles', maxCycles)
    p.set_param('timeStep', 0.1)
    p.set_param('spaceWidth', space_width) # i.e. 2cm x 2cm
    p.set_param('writeMediaLog', False)
    p.set_param('MediaLogRate', 1)
    p.set_param('writeFluxLog', False)
    p.set_param('FluxLogRate', 1)
    p.set_param("totalBiomassLogRate", 1)

    ## run the model, then save total biomass
    sim = c.comets(l, p)
    sim.working_dir = base

    try:
        sim.run()
    except:
        print(sim.run_output)

    biomass_data = sim.total_biomass
    biomass_data["growth_rate"] = growth_rate
    biomass_data = biomass_data.tail(2)

    liquid_sim_data_final = pd.concat([liquid_sim_data_final, biomass_data])
    print(i + 1, "of", max(expanded_grid.index.values) + 1)

liquid_sim_data_final.to_csv("/Users/abisesi/Desktop/PhD/Projects/Toxins-Growth-Rate/toxin-simulations/simulations_liquid/liquid_biomass_file_figure.csv")

import glob
for f in glob.glob("/Users/abisesi/Desktop/PhD/Projects/Toxins-Growth-Rate/toxin-simulations/simulations_liquid/totalbiomasslog*"):
    os.remove(f)

import glob
for f in glob.glob("/Users/abisesi/Desktop/PhD/Projects/Toxins-Growth-Rate/toxin-simulations/simulations_liquid/*.cmd"):
    os.remove(f)

import glob
for f in glob.glob("/Users/abisesi/Desktop/PhD/Projects/Toxins-Growth-Rate/toxin-simulations/simulations_liquid/*.txt"):
    os.remove(f)