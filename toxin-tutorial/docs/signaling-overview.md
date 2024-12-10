# An overview of the COMETS signaling capabilities

The COMETS platform is now equipped with two additional modeling capabilities: modeling single signal-reaction relationships, where a single metabolite changes the bounds on a single reaction, and modeling multiple signal-single reaction relationships, where several metabolites additively change the bounds on a single reaction. 

From these capabilities, a number of useful modeling scenarios can be generated:

* **Interspecies toxins / signaling**: One species produces a book-keeping metabolite that a second species can see via an exchange reaction and which signals changes in the second species. This is the basis of our manuscript.
* **Quorum-sensing**: A book-keeping metabolite is produced and secreted as part of the growth objective and serves as a signal that increases the bounds on a costly reaction.
* **Diauxie**: A nutrient metabolite influences the transport reaction of a second metabolite.
* **Multitoxins**: Multiple metabolites lower the upper bound on a species' biomass.

We discuss the basic structure of the `add_signal` and `add_multitoxin` functions which encode this ability. We also provide several examples of building signal-sensitive or signal-producing models and setting up well-mixed and spatial simulations appropriately. 

## The signaling function

There are two functions that allow users the ability to add signal production and response to their metabolic models: `add_signal`, which sets up a single metabolite-reaction relationship, and `add_multitoxin`, which sets up a multiple metabolite-single reaction relationship. The basic structure of the `add_signal` function looks like this in `cometspy`:

    model.add_signal(altered_reaction_id, metabolite_id, altered_bound, functional_relationship, [parms])

The biological significance of these parameters is:

* `altered_reaction_id`: The numeric ID of the metabolic reaction whose bounds will be altered by the uptake of the signaling compound. If modifying the maximum growth rate, the biomass growth objective should be used. When modeling -cidal toxins, the argument `'death'` can be used instead.
* `metabolite_id`: The index of the exchange reaction for the signaling metabolite. Models sense the external environment based on exchange reactions, which is why a reaction is specified rather than a metabolite. This value should therefore be numeric and correspond to the `EXCH_IND` of the appropriate metabolite.
* `altered_bound`: The reaction bound that is impacted by uptake of the signaling compound, either `'ub'` (upper bound) or `'lb'` (lower bound). Generally when modeling toxins that impact growth rate, `ub` should be used, although `lb` may be appropriate if users are interested in using a signaling compound to upregulate some metabolite process.
* `functional_relationship`: The signal-response curve detailing the shape of the functional relationship between the signaling compound concentration and the altered reaction bound. Three options are available: `'linear'`, `'bounded_linear'`, and `'generalized_logistic'`. `linear` and `bounded_linear` will model a fixed effect. The choice of functional relationship will determine the length and identities of the list of inputs to the final argument, `parms`.
* `parms`: The parameters for the functional relationship between signaling compound concentration and the reaction bound. The number of parameters and their meaning will vary based on the functional relationship chosen. The appropriate parameters should be included as a list, as indicated by the brackets.

## The signal-response curve

One essential decision when designing simulations with signaling compounds is the functional relationship between the concentration of the compound and the altered reaction. We provide the ability to model both linear and sigmoid relationships. The appropriate choice will depend strongly on the biological details of the system that the user is hoping to model, although there are several useful considerations, including the number of parameterizing arguments required to implement the different functional relationships.

If the functional relationship `linear` is chosen, two parameters are necessary for the `parms` argument, such that `parms = [slope, bound = 0]`:

* `slope`: The slope of the linear relationship between the concentration of the signaling compound and the bound of the altered reaction.
* `bound`: The maximum (or minimum, if using `lb`) amount of flux that the altered reaction (`altered_reaction_id`) can carry in the absence of toxin or signaling compound. When this parameter is not provided, it takes a default value of 0. 

If the functional relationship `bounded_linear` is chosen, four parameters are necessary for the `parms` argument, such that `parms = [baseline_value, conc_where_effect_starts, slope, conc_where_effect_saturates]`:

* `baseline_value`: The default response of the altered reaction at concentrations below the critical threshold `conc_where_effect_starts`. 
* `conc_where_effect_starts`: The signaling compound concentration at which the altered reaction begins to change. 
* `slope`: The slope of the linear relationship between the concentration of the signaling compound and the bound of the altered reaction. 
* `conc_where_effect_saturates`: The signaling compound concentration at which the reduction of the altered reaction saturates. 

If the functional relationship `generalized_logistic` is chosen, four parameters are necessary, and up to seven can be provided, for the `parms` argument, such that `parms = [left_asym, right_asym, growth_rate, starting_dose = 0, C = 1, Q = 1, v = 1]`:

* `left_asym`: The default response of the altered reaction at concentrations below the critical threshold, which in a dose-response logistic function is the left horizontal asymptote.
* `right_asym`: The signaling compound concentration at which the reduction of the altered reaction saturates, which in a dose-response logistic function is the right horizontal asymptote.
* `growth_rate`: The growth rate, or steepness, of the dose-response curve. 
* `starting_dose`: The signaling compound concentration at which the altered reaction begins to change, or move away from the left horizontal asymptote. If this parameter is not provided, it takes a default value of 0.

The parameters `C`, `Q` and `v` all take default value of 1, which should cover most use-cases of the generalized logistic function for signaling purposes. Unless users have a clear understanding of the dose-response curve they are interested in modeling and a strong rationale for setting these parameters themselves, it is advisable to only provide `left_asym`, `right_asym`, `growth_rate` and `starting_dose` parameters and use the default values for the remaining three. 

Note that the `bounded_linear` relationship is preferable over `linear` when users would like to prevent a reaction like biomass from becoming negative as a result of toxin. `generalized_logistic` will generally provide the most biologically-realistic relationship between model reactions and signaling compound concentrations, although there are many situations in which `bounded_linear` and `generalized_logistic` will not provide qualitatively divergent results. 

## Implementing multitoxins

Finally, users may chose a special case of signaling with `add_multitoxin`. The function adds a signaling relationship for multiple toxins, relying on a Hill-function reduction of (typically) the upper bound of a reaction. Users will find this function most useful when considering multiple, additively-functioning signals. The function arguments are very similar to the `add_signal` function, though multiple exchange IDs, along with compound-associated Km values and Hill coefficients, are required. 

    model.add_multitoxin(altered_reaction_id, [exch_ids], bound, vmax, [kms], [hills])

The parameters here are:

* `altered_reaction_id`: The numeric ID of the metabolic reaction whose bounds will be altered by the uptake of the signaling compounds.
* `exch_ids`: A list of exchange indexes for the signaling metabolite. 
* `bound`: The bound of the reaction to be altered, either `ub` or `lb`.
* `vmax`: The maximum bound of the reaction in the absence of signals.
* `kms`: A list of Km values (half-saturation concentrations) for the signals being incorporated.
* `hills`: A list of Hill coefficients for all the signals being incorporated.

## Useful notes for successful modeling

There are several important modeling limitations to consider when setting up signaling in a COMETS simulation.

* When users choose to model -cidal toxins with the `'death'` argument, instead of specifying whether the signal should change the upper or lower bound on the altered reaction, they can instead specify whether they want the toxin to be consumed (removed from the environment at a rate equal to cellular uptake) or not consumed, using `'met_unchanged'` or `'consume_met'`. 
* Users should take care when designing the signals to ensure that they are not designing impossible scenarios, such as the lower bound of a flux being greater than the upper bound. Plan to implement a variety of check to test that the upper bound of an altered reaction is not being reduced past the provided lower bound of that reaction and vice versa.  
* If users set up multiple signals to act on a single reaction without using the `add_multitoxin` function, the signal that was added last will override any signals previously added. The only exception to this is if two different signals affect different bounds on the same reaction.
* It is not advisable to set signals to affect exchange reactions. COMETS uses the bounds on exchange reactions for a variety of purposes, so changing their bounds is not a good idea. Instead, users can change *transport* reactions if they want to influence uptake using a signal.
* In order to be affected by a signal, a model must have an exchange reaction associated with that metabolite. It may need to be added to the model beforehand, but it is not necessary for the exchange reaction to be connected to any internal stoichiometry.
* COMETS uses Euler integration for signaled effects and death. This means simulation results will be sensitive to the timestep used. Users are advised in later portions of this tutorial to test how reducing the timestep influence the stability of the modeling results.