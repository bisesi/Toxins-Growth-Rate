# Building metabolic models with signalling compound production or response

The first step to running simulations with signalling compounds is building metabolic models that can accommodate the production of or response to these compounds. We will detail the model building process used in our paper, although this protocol can also be adapted to make changes to existing models rather than developing toy models. 

## Making toy models using `cobra`

The most straightforward way to generate toy metabolic models is to create a function that defines the desired models using a variety of functions and object types from `cobrapy`.

    def make_RPS_cobra_models_biomass_cost(growth_rate = 1., toxin_cost = 0.02, resistance_cost = 0.0, toxin_prod = 1.):
        # create metabolites for carbon and toxin
        carbon_e = Metabolite(id = "carbon_e",
                compartment = "e")
        carbon_c = Metabolite(id = "carbon_c",
                    compartment = "c")
        toxin_c = Metabolite(id = "toxin_c", compartment = "c")
        toxin_e = Metabolite(id = "toxin_e", compartment = "e")

        # create the reactions - biomass, carbon exchange and transport, and toxin exchange and transport
        #exchange reactions
        EX_carbon_e = Reaction(id = "EX_carbon_e",
                      lower_bound = -growth_rate, # growth rate
                      upper_bound = 1000.)
        EX_carbon_e.add_metabolites({carbon_e: -1})
        EX_toxin_e = Reaction(id = "EX_toxin_e",
                         lower_bound = -1000.,
                         upper_bound = 1000.)
        EX_toxin_e.add_metabolites({toxin_e: -1})

        # transport reactions
        carbon_transfer = Reaction(id = "carbon_transfer",
                          lower_bound = 0.,
                          upper_bound = 1000.)
        carbon_transfer.add_metabolites({carbon_e: -1,
                                carbon_c: 1})
        toxin_transfer = Reaction(id = "toxin_transfer",
                             lower_bound = -1000.,
                             upper_bound = 0.)
        toxin_transfer.add_metabolites({toxin_e: -1,
                                   toxin_c: 1})
        
        # biomass reactions
        Biomass = Reaction(id = "Biomass",
                  lower_bound = 0.,
                  upper_bound = 1000.)
        Biomass.add_metabolites({carbon_c: -1.})
        Biomass_producer = Reaction(id = "Biomass",
                  lower_bound = 0.,
                  upper_bound = 1000.)
        Biomass_producer.add_metabolites({carbon_c: -(1. + toxin_cost + resistance_cost),
                                     toxin_c: toxin_prod * (1. + toxin_cost + resistance_cost)})
        Biomass_resistant = Reaction(id = "Biomass",
                  lower_bound = 0.,
                  upper_bound = 1000.)
        Biomass_resistant.add_metabolites({carbon_c: -(1. + resistance_cost)})
        # generate the 3 model types with the appropriate reactions and objectives
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

This function takes four arguments (growth rate, the cost of toxin production, the cost of toxin resistance, and the amount of toxin produced per gram of biomass produced) to create three model genotypes (susceptibles, toxin producers, and resistant cheaters). Susceptible models have four metabolic reactions (carbon exchange, carbon transport, biomass production, and extracellular to intracellular toxin exchange), producers have five metabolic reactions (carbon exchange, carbon transport, biomass production, intracellular to extracellular toxin exchange, and toxin transport), and resistant cheaters have three metabolic reactions (carbon exchange, carbon transfer, and biomass production). 

When manipulating this function to generate toy models, particular attention should be paid to the construction of the biomass reactions, which will allow users to change the relative growth rates of each model (the default here is that they all have the same growth rate, aside from minimal penalities for the cost of toxin production or resistance) along with a variety of other assumptions. Here, toxin production is tied to producer biomass and carbon uptake:

```
Biomass_producer.add_metabolites({carbon_c: -(1. + toxin_cost + resistance_cost),
                                     toxin_c: toxin_prod * (1. + toxin_cost + resistance_cost)})
```

This setup defines toxin production in such a way that toxins are only created during growth and that the costs associated with toxin production (`toxin_cost`) or resistance to toxin (`resistance_cost`) vary directly with growth rate. It also establishes that toxin production is directly proportional to carbon uptake, regardless of growth rate, although the amount of toxin created per unit of carbon can be changed (`toxin_prod`). If desired, these constraints can be altered to fit the appropriate biological details of the signalling compounds being modeled. 

Notably, the function defined above also penalizes the biomass production of resistant cheaters to impose a cost of toxin resistance (`resistance_cost`), seen here:

```
Biomass_resistant.add_metabolites({carbon_c: -(1. + resistance_cost)})
```

Many of these modeling assumptions can be adjusted by changing the biomass reactions for any of the three models, or by altering the four arguments in the function call:

    growth_rate = 0.75 # growth rate for all three genotypes
    toxin_cost = 0.02 # 2% reduction in growth rate as a result of producing energetically expensive toxin
    resistance_cost = 0 # no reduction in growth rate as a result of being resistant to toxin
    toxin_prod = 15 # 15 mmol of toxin produced per gram of producer biomass

    producer, resistant, susceptible = make_RPS_cobra_models_biomass_cost(growth_rate = growth_rate, toxin_cost = toxin_cost, resistance_cost = resistance_ost, toxin_prod = toxin_prod)

## Making COMETS models from `cobra` models

Once `cobrapy` toy models have been constructed, they need to be converted into models that are recognizable by `cometspy`. Exactly how this is done will vary slightly for each type of model and its associated biological details.

    P = c.model(producer) # create the comets model from the cobra producer
    R = c.model(resistant)
    S = c.model(susceptible)

    P.obj_style = "MAX_OBJECTIVE_MIN_TOTAL" # set flux style as FBA (pFBA is also possible)
    R.obj_style = "MAX_OBJECTIVE_MIN_TOTAL"
    S.obj_style = "MAX_OBJECTIVE_MIN_TOTAL"

During the initial conversion, COMETS will erroneously assume that the biomass equation is an exchange for the susceptible and resistant models, because the reactions are unbalanced (unlike the producer, for which carbon uptake is balanced by toxin production). That will need to be fixed first. 

    # tell COMETS that the biomass reaction is not an exchange
    R.reactions.loc[R.reactions.REACTION_NAMES == "Biomass", "EXCH"] = False # set false
    index = R.reactions.loc[R.reactions.REACTION_NAMES == "Biomass", "EXCH_IND"].values[0]
    R.reactions.loc[R.reactions.REACTION_NAMES == "Biomass", "EXCH_IND"] = 0
    R.reactions.loc[R.reactions.EXCH_IND > index, "EXCH_IND"] -= 1

    S.reactions.loc[S.reactions.REACTION_NAMES == "Biomass", "EXCH"] = False
    index = S.reactions.loc[S.reactions.REACTION_NAMES == "Biomass", "EXCH_IND"].values[0]
    S.reactions.loc[S.reactions.REACTION_NAMES == "Biomass", "EXCH_IND"] = 0
    S.reactions.loc[S.reactions.EXCH_IND > index, "EXCH_IND"] -= 1

Next, we will need to add the signalling capability to the susceptible model, so that the uptake of toxin results in a proportionate reduction in growth rate. This can be done using a function from `cometspy` called `add_signal`. 

To add a signal, first retrieve the ID of the biomass and toxin exchange reactions, as these are important arguments for the `add_signal` function.

```
biomass_id = S.reactions.loc[S.reactions.REACTION_NAMES == "Biomass", "ID"].values[0]
toxin_exch_id = S.reactions.loc[S.reactions.REACTION_NAMES == "EX_toxin_e", "EXCH_IND"].values[0]
```

The `add_signal` function also requires a baseline growth rate in the absence of the signalling molecule, a threshold concentration of the signaling molecule where the effect on growth rate starts, the slope of the effect on growth rate, and the concentration of the signaling molecule where the effect on growth rate saturates. Adding a signal to the susceptible model will therefore look something like this:

```
gr_absense_of_toxin = growth_rate # baseline growth rate
toxin_conc_where_effect_starts = 0. # threshold concentration where effect starts
slope_of_toxin_effect = -growth_rate # slope of the relationship betwen toxin concentration and growth rate
toxin_conc_where_effect_saturates = 1. # threshold concentration where the effect ends
add_signal_parameter = 'bounded_linear'
S.add_signal(biomass_id, toxin_exch_id, 'ub', add_signal_parameter, 
             parms = [gr_absense_of_toxin, 
                      toxin_conc_where_effect_starts, 
                      slope_of_toxin_effect,
                     toxin_conc_where_effect_saturates])
```
Using these arguments will create a signal such that the growth rate of the susceptible model will reduce in proportion to the concentration of toxin in the environment, with the effect on growth rate occurring whenever toxin is present and saturating when the toxin concentration reaches 1. The `add_signal_parameter` tunes the shape of relationship between toxin concentration and biomass. Here, the `"bounded_linear"` argument ensures that the upper bound of growth rate does not become negative. 

## Important notes on toxins and signalling 

While this tutorial provides an example of toxin production, the `add_signal` function can accommodate a variety of approaches for modeling all types of signalling compounds. Manipulating the arguments of the function will determine the metabolic consequences of a signalling molecule for a focal model. The function has eight main arguments:

    model.add_signal(altered_reaction_id, metabolite_id, altered_bound, functional_relationship, parms = [baseline_value, conc_where_effect_starts, slope_of_effect, conc_where_effect_ends])

The biological significance of these parameters is:

* `altered_reaction_id` = the numeric ID of the metabolic reaction whose bounds will be altered by the uptake of the signalling compound, generally the biomass reaction
* `metabolite_id` = the ID of the signalling compound exchange reaction
* `altered_bound` = the bound that is impacted by signalling compound uptake, either `'ub'` (upper bound) or `'lb'` (lower bound)
* `functional_relationship` = the shape of the functional relationship between the signalling compound concentration and the altered reaction bound, with three options available: `'linear'`, `'bounded_linear'`, and `'generalized_logistic'`
* `baseline_value` = the default value of the flux through the altered reaction
* `conc_where_effect_starts` = the concentration of the signalling compound where the effect on the reaction starts
* `slope_of_effect` = the slope of the functional relationship between the signalling compound concentration and the altered reaction bound
* `conc_where_effect_ends` = the concentration of the signalling compound where the effect on the reaction saturates

