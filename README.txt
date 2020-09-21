--- GAFFER ---

The gaffer library includes a genetic algorithm meant to partition time series data into spatial clusters depending on relationships with covariates.

For example, suppose we have weekly counts of human malaria cases in a number of districts, and we suspect that malaria cases in some districts respond strongly to recent precipitation, while risk in other districts does not respond much at all to rain. We would like to put these rain-sensitive districts into one cluster and the less-sensitive districts into another.

At the same time, we are also interested in determining which covariates are informative of risk - is precipitation even predictive? Or perhaps we should be looking at soil moisture? Or both?

In what follows, we describe what is necessary to run gaffer and what outputs the users will obtain throughout the run.

--- BEFORE THE RUN ---

Most data are organized on the basis of districts and dates. Districts may be countries, counties, zip codes, etc. as long as they are disjoint spatial units. Dates indicate the beginning of the period for which the dependent variable is reported; e.g. the date beginning the week for which cases are reported.

Data will include:

- One dependent variable (assumed to be counts of cases of malaria or some other disease).
- Multiple independent variables (assumed to be environmental variables measured at the same times and locations as the disease data).
- Any other variables that are referenced in the model formula (e.g. )
- Each data record references a given location (district, county, etc.) at a given date.
- Best results when we have a panel data format, however, the algorithm can handle missing data values.
- User determines the geographic area prior to the run, and provides appropriate data for that area.



The user provides the following data arguments:
Data frame of case counts - has space/time fields.
Data frame of environmental variables - has space/time fields.
Spatial data - adjacency matrix (R list object)

The user passes the model formula (including the model family) as a parameter.

The user can either accept default parameters for the GA run or specify their own parameters.
Number of generations
Number of initial clusters
Number of covariates
Individuals per generation
Mutation probabilities
Whether or not there is an alpha selection
Model measure (we’re using only AIC now, but this is something that we would like to modify in the future to provide more options)

The user can specify two options for variable selection
The user provides a set of predictor variables in the input file and specifies a maximum number of environmental variables (m) that can be included in the final model. The algorithm will explore models with different combinations of environmental variables and select the best model containing m or fewer variables.
The users specifies the number and identify of the predictor variables to be used in the model (e.g., 3 variables: lst_day, prec, and ndvi). The variables are used in all models, such that the only differences among the models are the assignments of clusters.

--- DURING THE RUN ---

The user has to wait a while to get results. What information do we want to provide them during the model run?

The user can monitor the progress of the run. For example, if the model does not seem to be converging they may want to abort the run and try something different.

On the screen
Generation number
Model measure
???
Other possibilities
Continually update the .csv output files - write to them every generation or every other generation. Append to the existing files, make multiple copies. But this is very inefficient. If we implement this, probably don’t make it the default.
Generate a convergence graph every iteration (or every n iterations) and write to the base graphics device. Better to do this with base graphics rather than ggplot.

The user can bail out of the run with Esc or Ctrl-C if they need to. It would be nice to incrementally save the output files after every iteration in case the user bails out and wants to go back and examine them more closely.

--- AFTER THE RUN ---

The user wants to be able to determine whether the algorithm has converged. To do that, the user will want to graph changes in the model measure (and possibly other variables) as a function of generation. The user may also want to be able to compare the spatial details of the “best” model with other competing models (top 10?).

Output tables saved as csv.

List of models per generation with
Seeds
Generation
Mutation
Model number
Model measure
Selection probability
Covariates
Cluster seeds

Details of the best model (Ram - provide more details here)

Cluster assignments per woreda for the best model (Ram provide more details here)


Additional features:

The user should have the capability to
Provide a specific starting configuration for the GA
Supply the details of a previous “best” model as the starting configuration for the GA.



