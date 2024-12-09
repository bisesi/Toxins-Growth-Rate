# Running simulations

Once metabolic models with signaling compound production and/or response have been built, it's time to run some simulations.

## Well-mixed environments

Toxin production can be simulated in a well-mixed environment using our three metabolite models: susceptibles, producers, and resistant cheaters. Recall that if producer models are not relevant for the simulations being run, signaling compounds instead should be added into the starting media at a constant or variable level, as appropriate.

    S.initial_pop = [0,0,1e-10] # initialize populations with the desired biomass, no x and y coordinates
    R.initial_pop = [0,0,1e-10] # resistant
    P.initial_pop = [0,0,1e-10] # producer

    l = cometspy.layout([P, R, S]) # set up the COMETS layout
    l.set_specific_metabolite("carbon_e", 4e-7) # add the carbon source into the layout

    p = cometspy.params() # specify the simulation parameters
    p.set_param("defaultVmax", 1.) # set Vmax
    p.set_param("defaultKm", 0.000001) # set Km
    p.set_param('maxCycles', 300) # set the total number of simulation hours
    p.set_param('timeStep', 0.1) # set the timestep of each simulation
    p.set_param('spaceWidth', 0.02) # set the size of the 2D environment, i.e. 2cm x 2cm
    p.set_param('writeMediaLog', True) # if True, COMETS will provide the log of metabolite concentrations over time
    p.set_param('MediaLogRate', 1) # timestep at which to record changes in metabolite concentrations
    p.set_param('writeFluxLog', True) # if True, COMETS will provide the log of metabolic fluxes for each model over time
    p.set_param('FluxLogRate', 1) # timestep at which to record changes in model flux
    p.set_param("totalBiomassLogRate", 1) # timestep at which to record changes in biomass

Once the starting population and simulation parameters have been set to the user's desires, simulations can be run.

```
sim = cometspy.comets(l, p)
sim.run()
```

Users can view the results of liquid simulations by examining the dataframes `sim.total_biomass` and `sim.media`. Note that the maximum number of cycles provided here may not be sufficient for the simulation to reach stationary phase depending on the growth rate provided and other parameters. Users should look carefully at the results of their simulations to ensure that they have gathered the information intended. 

In well-mixed simulations, five parameters are generally most consequential for the outcome of competition: the growth rate cost of toxin production, the toxin production rate, the total amount of carbon provided in the environment, the total size and relative genotype frequencies of the initial population, and the parameter `spaceWidth`. The `spaceWidth` parameter influences the concentration of the toxin and is therefore an important value when modeling signalling compounds. Smaller `spaceWidth` values will result in higher effective concentrations of toxins or other signalling molecules. 

## Spatially-structured environments

Signaling compound production can also be simulated in 2D space using the same three metabolite models: susceptibles, producers, and resistant cheaters. Again, if compound production is not relevant for simulations, the compound will need to be added into the media in some other way. 

### Setting up the initial location of colonies

One of the first and most important decisions that arises when running spatial simulations is how to initialize the starting locations of colonies. We detail the strategy taken in our paper, but note that COMETS has additional capabilities to set locations deterministically or randomly. 

To randomly generate starting locations (x,y coordinates) for colonies, users can define the following function:

    def pick_unique_locations(width, height, n, edge_space = 0):
        locs = []
        while len(locs) < n:
            loc = (random.randrange(edge_space, width - edge_space),
                random.randrange(edge_space, height - edge_space))
            if loc not in locs:
                locs.append(loc)
        return(locs)

The use of the `random` package from `Python` also necessitates that users set a `spatial_seed` for reproducibility:

```
random.seed(spatial_seed) 
```

This will ensure that the `pick_unique_locations` function uses the same initial starting locations each time the same `spatial_seed` is set. 

Once the function is defined and the `spatial_seed` for the simulation has been assigned, unique starting locations for any desired number of colonies can be generated. The relative ratios of the three genotypes will also be determined by the number of starting locations assigned to each genotype. 

    grid_size = [50,50] # size of the 2D grid from which locations should be picked
    number_of_colonies = 30 # number of initial colonies for which x,y coordinates are needed
    edge_space = 3 # distance from the grid border, beyond which coordinates should not be picked
    locs = pick_unique_locations(grid_size[0], grid_size[1], number_of_colonies, edge_space)

    producer_locs = locs[0:5] # the first four coordinates are assigned to be producers
    resistant_locs = locs[5:10] # the next four coordinates are assigned to be resistants
    susceptible_locs = locs[10:] # the remaining coordinates are assigned to be susceptibles

    initial_biomass = 1.e-10 # biomass for each colony
    S.initial_pop = [[loc[0], loc[1], initial_biomass] for loc in susceptible_locs] # assign the susceptible model to all susceptible colony locations with the appropriate starting biomass
    P.initial_pop = [[loc[0], loc[1], initial_biomass] for loc in producer_locs]
    R.initial_pop = [[loc[0], loc[1], initial_biomass] for loc in resistant_locs]  

    l = c.layout([P, R, S]) # set up the initial COMETS layout

It is generally desirable to save the starting locations of colonies in order to do any spatially-explicity analyses later. There are many approaches to accomplish this. For example, the strain genotype and x,y coordinates can be assigned to a dataframe and exported as a csv. Users should also consider saving the `spatial_seed` value set for the generation of starting locations in order to make their analyses reproducible. 

### Tracking the biomass of individual colonies

One common challenge of spatially-structured simulations is that the default settings do not accommodate tracking the biomass of individual colonies. Instead, simulations can track the total biomass of each model type (susceptible, producer, resistant cheater) and the biomass of each model type at each x,y coordinate. Because biomass diffuses, many x,y locations in a grid will have a non-zero biomass value by the end of a simulation, and while the type of model that that biomass corresponds to can be discerned, which initial colony the biomass should be counted toward is not possible to ascertain in most cases. 

There are a few approaches to circumvent this limitation. The most straightforward involves assigning a unique metabolic model to each location, so that the total biomass file reflects the growth of every individual colony. This can be accomplished by making multiple individual COMETS models for all genotypes or models of interest.

    grid_size = [50,50]
    number_of_colonies = 30 
    edge_space = 3 
    locs = pick_unique_locations(grid_size[0], grid_size[1], number_of_colonies, edge_space)
    producer_locs = locs[0:5]
    resistant_locs = locs[5:10]
    initial_biomass = 1e-10 # total initial biomass for each colony

    producer_models = []
    for l in range(0, len(producer_locs)):
        P = c.model(producer) # take the cobra model and convert it into a comets model

        P.obj_style = "MAX_OBJECTIVE_MIN_TOTAL" # set the fba objective 

        loc = producer_locs[l] # get the starting location
        P.initial_pop = [loc[0], loc[1], initial_biomass] # assign starting location and biomass
        producer_models.append(P) # add to list of other producer models with different locations

    resistant_models = []
    for k in range(0, len(resistant_locs)):
        R = c.model(resistant)

        R.obj_style = "MAX_OBJECTIVE_MIN_TOTAL"

        # fix COMETS thinking biomass reactions with only reactants are exchanges
        R.reactions.loc[R.reactions.REACTION_NAMES == "Biomass", "EXCH"] = False
        index = R.reactions.loc[R.reactions.REACTION_NAMES == "Biomass", "EXCH_IND"].values[0]
        R.reactions.loc[R.reactions.REACTION_NAMES == "Biomass", "EXCH_IND"] = 0
        R.reactions.loc[R.reactions.EXCH_IND > index, "EXCH_IND"] -= 1

        loc = resistant_locs[k]
        R.initial_pop = [loc[0], loc[1], initial_biomass]
        resistant_models.append(R)

This will generate a list of `cometspy` producer or resistant models with the associated x,y coordinate and initial biomass for each colony. These models can then be concatenated into a longer list prior to the start of simulations in order to create the desired layout. 

```
 l = c.layout(producer_models + resistant_models)
```

When simulations are initialized in this way, the total biomass file will contain each colony as a column (i.e. producer1, producer2), with each row as the biomass of that colony at a given timestep. This type of analysis makes it especially important that users save the starting locations of colonies. 

### Running the simulations

Once users have initialized the starting locations of colonies, it's time to run simulations. The process is similar to that for well-mixed environments, with the addition of a handful of important variables.

    space_width = 0.01 # simulation environment size i.e. 1 cm x 1 cmm
    max_cycles = 30 # initial guess for the time to stationary, should be adjusted as needed
    dt = 0.25 * ((space_width)**2) / metabolite_diff # dt should be determined based on diffusion rates and environment size
    metabolite_diff = 5e-6 # default rate of metabolite (carbon) diffusion
    toxin_diff = 5e-7 # default rate of toxin diffusion
    grid_size = [50,50] # coordinate size of the spatial grid to simulate

    p = c.params() # set up simulation parameters
    p.set_param("defaultVmax", 1.)
    p.set_param("defaultKm", 0.000001) 
    p.set_param('maxCycles', max_cycles)
    p.set_param('timeStep', dt) # this needs to be calculated relative to diffusion and spatial simulation size
    p.set_param('spaceWidth', space_width) # this helps determine timeStep
    p.set_param('writeMediaLog', True)
    p.set_param('writeBiomassLog', False) # if True, provides a log of biomass at each x,y location in the environment grid
    p.set_param('minSpaceBiomass', 1.e-15)
    p.set_param('defaultDiffConst', metabolite_diff) # rate of metabolite diffusion
    p.set_param('flowDiffRate', 3.e-9)
    p.set_param('writeFluxLog', False) 
    p.set_param("totalBiomassLogRate", 1)
    p.set_param("MediaLogRate", 1)
    p.set_param("numDiffPerStep", 1) # rate of diffusion per time step

    l = c.layout(producer_models + resistant_models) # make the layout with producer and resistant models
    l.grid = grid_size # set up the grid based on the appropriate size
    l.set_specific_metabolite("carbon_e", 4.e-7) # initialize the starting amount of carbon at each x,y coordinate
    l.media.loc[l.media.metabolite == "toxin_e", "diff_c"] = toxin_diff # add the rate of toxin diffusion

    # run the simulation
    sim = c.comets(l, p)
    sim.run(delete_files = False)

Users can examine the outputs of `sim.total_biomass` and `sim.media` to evaluate their simulations. It is almost certain that individual adjustments will need to be made to ensure that colonies reach stationary or that the rate of diffusion is appropriate depending on the other desired parameters like growth rate, etc. 

### Troubleshooting appropriate simulation durations

For successful spatial simulations, users should be attentive to the timestep (`timeStep` or `dt`), the rate of diffusion (`defaultDiffConst` and `numDiffPerStep`), and total simulation duration (`maxCycles`). 

To test that the timestep `dt` (`p.set_param('timeStep', dt)`) is appropriate, users should consult the `sim.run_ouput` following an initial short simulation. `run_ouput` will provide information about the diffusion constants, in particular raising a warning if the diffusion constants are unstable. The file will also suggest the appropriate timestep to fix the problem. By calculating `dt` with the `spaceWidth` and `metabolite_diff` variables, users should be able to largely offset this issue. However, it is worth verifying before running computationally expensive simulations with a higher number of `maxCycles` in order to make sure that estimations are as robust as possible. 

Additionally, the appropriate number of `maxCycles` required to reach the cessation of growth (i.e. stationary phase) is likely to be highly dependent on simulation-specific parameters like diffusion constants and growth rate. This should be adjusted as needed and may require several rounds of running simulations and checking the final timesteps to see if biomass continues to change. 