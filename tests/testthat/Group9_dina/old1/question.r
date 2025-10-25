setCDM <- function(
    model, population_distribution = NULL, q_matrix = NULL, prior_pi = NULL,
    rule = c("DINA", "DINO")) {
    if (missing(rule)) {
        warning(
            "`rule` not supplied; defaulting to 'DINA'. ",
            "Use rule = 'DINO' if needed."
        )
        rule <- "DINA"
    } else {
        rule <- match.arg(rule)
    }

    new("cdm",
        model = model,
        population_distribution = population_distribution,
        q_matrix = q_matrix,
        prior_pi = prior_pi,
        rule = rule
    )
}
This above function sets up an instance for the cognitive diagnostic model. Downstream of the process, the model optimisation will need q_matrix and prior_pi. How can I have 
this function to set up a default q_matrix and prior_pi? Specifically, the 'model'
instance the function got is something like: 

model <- BuildModel(
    p_map = list(
        guess1 = "1", guess2 = "1", guess3 = "1",
        slip1 = "1", slip2 = "1", slip3 = "1"
    ),
    factors = NULL,
    constants = NULL,
    match_map = NULL,
    accumulators = NULL,
    type = "cdm",
    verbose = TRUE
)
It has the following slots:
slotNames(model)
 [1] "parameter_map"               "accumulators"               
 [3] "factors"                     "match_map"                  
 [5] "constants"                   "cell_names"                 
 [7] "parameter_x_condition_names" "model_boolean"              
 [9] "pnames"                      "npar
 
One can infer the number of items by looking up the slot, model@parameter_map:

$guess1
[1] "1"

$guess2
[1] "1"

$guess3
[1] "1"

$slip1
[1] "1"

$slip2
[1] "1"

$slip3
[1] "1",

which indicates that three items were tested (assumed), because the CDM has two 
types of parameters: guess and slip, and looking at the list, one can know that three guesses and three slips were pressumed, so the number of items is 3.


Although we will still need the user to tell us the number of pressumed skills (or sometimes called attributes), we may presume the smallest, reasonable number. That is, 2. (Presuming 1 would be meaningless to do the modelling.) And we can set up a default Q matrix like:

Q <- matrix(c(
    1, 0,
    0, 1,
    1, 1
), ncol = 2, byrow = TRUE)
colnames(Q) <- c("A1", "A2")

n_item <- nrow(Q)
n_skill <- ncol(Q)
n_profile <- 2^(n_skill)
pi_uniform <- rep(1 / n_profile, n_profile)

So does the default prior_pi (named as pi_uniform in the above code.)

Thus, my question is that how can I adjust the function
setCDM <- function(
    model, population_distribution = NULL, q_matrix = NULL, prior_pi = NULL,
    rule = c("DINA", "DINO")) {
    if (missing(rule)) {
        warning(
            "`rule` not supplied; defaulting to 'DINA'. ",
            "Use rule = 'DINO' if needed."
        )
        rule <- "DINA"
    } else {
        rule <- match.arg(rule)
    }

    new("cdm",
        model = model,
        population_distribution = population_distribution,
        q_matrix = q_matrix,
        prior_pi = prior_pi,
        rule = rule
    )
}
, to do the job described above and also give a message to the user to let them know that if they want they can enter a user defined q_matrix, prior_pi and rule themselves. These should not be warning because I do not want to let the user thinks that these messages will make the modelling wrong or does not work.