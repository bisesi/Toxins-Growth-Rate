# Building metabolic models with signaling compound production or response

The first step to running simulations with signaling compounds is building metabolic models that can accommodate the production of or response to these compounds. We will detail several model-building processes, including the approach taken in our paper.

## Adding new signals to existing models

Many users will be most interested in adding the signaling capability to existing genome-scale metabolic models. We will walk through this process by using an existing model as an example.

    E_cobra = cobra.io.load_model('textbook') # basic IJO1366 model

As with any COMETS simulation, `cobra` models need to be converted to COMETS models prior to usage. This is a straightforward process in the case of existing models.

    E = c.model(E_cobra)
    E.open_exchanges() # ensures that models can take up any nutrients provided in the COMETS media
    E.obj_style = "MAX_OBJECTIVE_MIN_TOTAL" # type of flux solution

Setting up a COMETS simulation also involves specifying the initial populations, creating a layout, and seeding the media with the relevant components. More information about this process is available in the COMETS documentation, as many of the specific details will vary significantly based on the desires of the user.

    E.initial_pop = [0, 0, 1.e-7] # the first two values are x, y coordinates, third is gDW
    layout = c.layout(E)
    E_media = {'co2_e': 1000.0,
             'h_e': 1000.0,
             'h2o_e': 1000.0,
             'nh4_e': 1000.0,
             'o2_e': 1000.0,
             'pi_e': 1000.0}
    for key, value in E_media.items():
        layout.set_specific_metabolite(key, value)
    layout.set_specific_metabolite("glc__D_e", 5e-4)

Once users have decided on many of the core parameters for their simulations, they may want to add a signaling compound. In this example, we will consider an undegradable toxin whose purpose is to reduce the growth of *E. coli* in a dose-dependent way. This toxin will exist in the media from the outset of simulations. 

To make a toxin, we will assign a signal to the *E. coli* model that alters a reaction in that model based upon the concentration of an environmental metabolite. We will make a new metabolite for this purpose, although this isn't strictly necessary, as we will see later. 

When setting a new metabolite as a toxin or signaling compound, the affected model will need to be altered to have an exchange reaction for that metabolite. We will use functions from `cobra` to do that.

    from cobra import Metabolite, Reaction
    toxin_e = Metabolite(id = "toxin_e", compartment = "e") # create a toxin metabolite in the extracellular environment
    EX_toxin_e = Reaction(id = "EX_toxin_e", lower_bound = -1000, upper_bound = 1000) # create an exchange reaction for the toxin 
    EX_toxin_e.add_metabolites({toxin_e: -1}) # add toxin metabolite and stoichiometric coefficient, which in this case is -1
    E_cobra.add_reactions([EX_toxin_e]) # add the exchange reaction to E coli

Now, the *E. coli* metabolic model can take up the toxin metabolite, which is essential for its sensitivity to the compound. However, we have yet to define how the compound will impact the growth rate or survival of *E. coli*, which we will do with the `add_signal` function. This is the essential innovation of this paper, and the arguments of the function are discussed in greater detail later in the tutorial. Here, we are interested in defining a toxin that will impact the growth rate of *E. coli* by reducing the upper bound on biomass without allowing biomass production to become negative. These changes need to be made to the COMETS version of the model.

    response = E.reactions.loc[E.reactions.REACTION_NAMES == "Biomass_Ecoli_core",:].ID.item() # reaction altered by the toxin
    signal_exch = E.reactions.loc[E.reactions.REACTION_NAMES == "EX_toxin_e",:].EXCH_IND.item() # exch_id of the toxin
    E.add_signal(response, signal_exch, 'ub', 'bounded_linear', parms = [1,0.2,-0.2,5.2])

The biomass of *E. coli* will now be reduced linearly based on the extracellular concentration of the toxin, with the effect starting at a concentration of 0.2mmol and saturating (i.e. reducing *E. coli*'s growth rate to 0) at a concentration of 5.2 mmol, with a slope of -0.2. These specific parameters can be changed to model different functional relationships, and more details about these values - and how to choose them - are discussed later.

We are now prepared to run a simulation with an *E. coli* model that will be sensitive to a toxin. In this case, we will need to add that toxin into the extracellular environment as an additional media component prior to running our simulation.

    toxin_mmol = 0.6 # set the concentration of toxin in the environment
    layout.set_specific_metabolite('toxin_e', toxin_mmol) # add toxin into the environment

## Turning existing metabolites into signaling compounds

In some cases, instead of creating an entirely new metabolite that serves as a signaling compound, users will want models to respond to existing metabolites. We detail that process here, based on the textbook *E. coli* model.

    model = cobra.io.load_model('textbook')

We can consider a situation where the extracellular build-up of the existing metabolite acetaldehyde (`EX_acald_e`) reduces the function of phosphofructokinase (`PFK`). The basics of the model set up will not change much from our previous example.

    # create COMETS model
    m = c.model(model)
    m.open_exchanges()
    m.initial_pop = [[0, 0, 0.01]]  

    # add signaling relationship without adding a new metabolite
    PFK_num = m.reactions.loc[m.reactions.REACTION_NAMES == "PFK", "ID"].values[0]
    acald_exch_id = m.reactions.loc[m.reactions.REACTION_NAMES == "EX_acald_e", "EXCH_IND"].values[0]

In this case, we will need to be more sensitive to the maximum flux through the `PFK` reaction in order to appropriately reduce it through the build up of acetaldehyde. To find the maximum flux, we can solve using the `cobra` model:

    sol = model.optimize()
    max_ub = sol.fluxes["PFK"]
    print(max_ub)
    intercept = max_ub

This value will be useful for creating our final signal. Because we are reducing the function of a metabolic reaction rather than growth rate, we will also change the linear relationship; the flux through `PFK` can be negative, so we will use the `linear` function_relationship rather than `bounded_linear`:

    slope = -1000 # this is a fairly arbitrary relationship for the slope between acetaldehyde concentration and flux through PFK
    m.add_signal(PFK_num, acald_exch_id, "ub", "linear", [slope, max_ub])

Here, we have defined a signal such that the maximum flux through PFK `max_ub` will reduce linearly with a slope of -1000 as the concentration of acetaldehyde increases. We can run our simulations in a couple of ways. First, as before, we will observe the consequences of a constant concentration of extracellular signaling compound.

    l = c.layout(m)
    l.set_specific_metabolite("glc__D_e", 5.)
    l.set_specific_metabolite("nh4_e", 1000.)
    l.set_specific_metabolite("pi_e", 1000.)
    l.set_specific_metabolite("acald_e", 0.25)

Using this approach, our signaling compound is present in the media at an initial concentration that will be reduced over simulation time as it is taken up by *E. coli*. 

We may also be interested in a signaling compound that is added dynamically, whether to replicate pulses by a producer or other scenario. The simplest way to do this is by using several available tools in COMETS, one of which is the ability to "refresh" compounds, adding a certain mmol of compound per spatial box per hour. 

    l.set_specific_metabolite("acald_e", 0.) # reset
    l.set_specific_refresh("acald_e", 0.25)

## Making toy producer models using `cobra`

Finally, users may want to make toy models, such as the ones used in our publication. Functionally, these models will not significantly differ from existing models. The basic premise of their construction is identical to the previous examples. Their benefit is that they are more easily manipulated, since every exchange and reaction is chosen and specifically defined by the user, and they are fractionally less computationally expensive due to having only a handful of fluxes. 

The most straightforward way to generate toy metabolic models is to create a function that defines the desired models using a variety of functions and object types from `cobrapy`. So far, we have only seen how to add signaling compounds that exist in the extracellular environment and impact a strain growing in isolation. This function provides a way to create a metabolic model that will produce the signaling compound. This same basic process can be applied to existing models, if users would like to add toxin production to existing models. 

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

When manipulating this function to generate toy models, particular attention should be paid to the construction of the biomass reactions, which will allow users to change the relative growth rates of each model (the default here is that they all have the same growth rate, aside from minimal penalities for the cost of toxin production or resistance) along with a variety of other assumptions. Here, we have tied toxin production to producer biomass and carbon uptake, which differs from previous examples, where signaling compounds already existed in the media:

```
Biomass_producer.add_metabolites({carbon_c: -(1. + toxin_cost + resistance_cost),
                                     toxin_c: toxin_prod * (1. + toxin_cost + resistance_cost)})
```

This setup defines toxin production in such a way that toxins are only created during growth and that the costs associated with toxin production (`toxin_cost`) or resistance to toxin (`resistance_cost`) vary directly with growth rate. It also establishes that toxin production is directly proportional to carbon uptake, regardless of growth rate, although the amount of toxin created per unit of carbon can be changed (`toxin_prod`). If desired, users can alter these constraints to fit the appropriate biological details of the signaling compounds being modeled. 

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

## Making COMETS models from `cobra` toy models

As we have already shown, once `cobrapy` toy models have been constructed, they need to be converted into models that are recognizable by `cometspy`. Exactly how this is done will vary slightly for each type of model and its associated biological details, and is somewhat more complicated when using toy models than pre-existing ones.

    P = c.model(producer) # create the comets model from the cobra producer
    R = c.model(resistant)
    S = c.model(susceptible)

    P.obj_style = "MAX_OBJECTIVE_MIN_TOTAL" # set flux style as FBA (pFBA is also possible)
    R.obj_style = "MAX_OBJECTIVE_MIN_TOTAL"
    S.obj_style = "MAX_OBJECTIVE_MIN_TOTAL"

Diverging from curated models, when using a toy model, during the initial conversion, COMETS will erroneously assume that the biomass equation is an exchange for the susceptible and resistant models, because the reactions are unbalanced (unlike the producer, for which carbon uptake is balanced by toxin production). That will need to be fixed first. 

    # tell COMETS that the biomass reaction is not an exchange
    R.reactions.loc[R.reactions.REACTION_NAMES == "Biomass", "EXCH"] = False # set false
    index = R.reactions.loc[R.reactions.REACTION_NAMES == "Biomass", "EXCH_IND"].values[0]
    R.reactions.loc[R.reactions.REACTION_NAMES == "Biomass", "EXCH_IND"] = 0
    R.reactions.loc[R.reactions.EXCH_IND > index, "EXCH_IND"] -= 1

    S.reactions.loc[S.reactions.REACTION_NAMES == "Biomass", "EXCH"] = False
    index = S.reactions.loc[S.reactions.REACTION_NAMES == "Biomass", "EXCH_IND"].values[0]
    S.reactions.loc[S.reactions.REACTION_NAMES == "Biomass", "EXCH_IND"] = 0
    S.reactions.loc[S.reactions.EXCH_IND > index, "EXCH_IND"] -= 1

Next, we will need to add the signaling capability to the susceptible model, so that the uptake of toxin results in a proportionate reduction in growth rate. 

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
Using these arguments will create a signal such that the growth rate of the susceptible model will reduce linearly in proportion to the concentration of toxin in the environment, with the effect on growth rate occurring whenever toxin is present and saturating when the toxin concentration reaches 1. Here, as we have previously seen, the `"bounded_linear"' argument ensures that the upper bound of growth rate does not become negative. 

## Details on the `add_signal` function

While this tutorial provides several examples of toxic metabolites, the `add_signal` function can accommodate a variety of approaches for modeling all types of signaling compounds. Manipulating the arguments of the function will determine the metabolic consequences of a signaling molecule for a focal model. The function has five main arguments:

    model.add_signal(altered_reaction_id, metabolite_id, altered_bound, functional_relationship, [parms])

The biological significance of these parameters is:

* `altered_reaction_id`: The numeric ID of the metabolic reaction whose bounds will be altered by the uptake of the signaling compound. If modifying the maximum growth rate, the biomass growth objective should be used. When modeling -cidal toxins, the argument `'death'` can be used instead.
* `metabolite_id`: The index of the exchange reaction for the signaling metabolite. Models sense the external environment based on exchange reactions, which is why a reaction is specified rather than a metabolite. This value should therefore be numeric and correspond to the `EXCH_IND` of the appropriate metabolite.
* `altered_bound`: The reaction bound that is impacted by uptake of the signaling compound, either `'ub'` (upper bound) or `'lb'` (lower bound). Generally when modeling toxins that impact growth rate, `ub` should be used, although `lb` may be appropriate if users are interested in using a signaling compound to upregulate some metabolite process.
* `functional_relationship`: The signal-response curve detailing the shape of the functional relationship between the signaling compound concentration and the altered reaction bound. Three options are available: `'linear'`, `'bounded_linear'`, and `'generalized_logistic'`. `linear` and `bounded_linear` will model a fixed effect. The choice of functional relationship will determine the list of inputs to the final argument, `parms`.
* `parms`: The parameters for the functional relationship between signaling compound concentration and the reaction bound. The number of parameters and their meaning will vary based on the functional relationship chosen. The appropriate parameters should be included as a list, as indicated by the brackets.

If the functional relationship `linear` is chosen, two parameters are necessary for the `parms` argument, such that `parms = [slope, bound]`:

* `slope`: The slope of the linear relationship between the concentration of the signaling compound and the bound of the altered reaction.
* `bound`: The maximum (or minimum, if using `lb`) amount of flux that the altered reaction (`altered_reaction_id`) can carry in the absence of toxin or signaling compound. 

If the functional relationship `bounded_linear` is chosen, four parameters are necessary for the `parms` argument, such that `parms = [baseline_value, conc_where_effect_starts, slope, conc_where_effect_saturates]`. The same basic formulation applies for `generalized_logistic` as well:

* `baseline_value`: The default response of the altered reaction at concentrations below the critical threshold `conc_where_effect_starts`. 
* `conc_where_effect_starts`: The signaling compound concentration at which the altered reaction begins to change. 
* `slope`: The slope of the linear relationship between the concentration of the signaling compound and the bound of the altered reaction. For `generalized_logistic` the slope of the curve between the upper and lower response limits of the dose-response curve.
* `conc_where_effect_saturates`: The signaling compound concentration at which the reduction of the altered reaction saturates. 

Note that the `bounded_linear` relationship is preferable over `linear` when users would like to prevent a reaction like biomass from becoming negative as a result of toxin. `generalized_logistic` will generally provide the most biologically-realistic relationship between model reactions and signaling compound concentrations, although there are many situations in which `bounded_linear` and `generalized_logistic` will not give qualitatively divergent results. 

Finally, users may chose a special case of signaling with `add_multitoxin`. The function adds a signaling relationship for multiple toxins, relying on a Hill-function reduction of (typically) the upper bound of a reaction. Users will find this function most useful when considering multiple, additively-functioning signals. The function arguments are very similar to the `add_signal` function, though multiple exchange IDs, along with Km values and Hill coefficients, are required. 

    model.add_multioxin(altered_reaction_id, [exch_ids], bound, vmax, [kms], [hills])

The parameters here are:

* `altered_reaction_id`: The numeric ID of the metabolic reaction whose bounds will be altered by the uptake of the signaling compounds.
* `exch_ids`: A list of exchange indexes for the signaling metabolite. 
* `bound`: The bound of the reaction to be altered, either `ub` or `lb`.
* `vmax`: The maximum bound of the reaction in the absence of signals.
* `kms`: A list of Km values (half-saturation concentrations) for the signals being incorporated.
* `hills`: A list of Hill coefficients for all the signals being incorporated.

