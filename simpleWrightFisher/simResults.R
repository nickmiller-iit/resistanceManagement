# S3 class and methods for working with simulations from SLiM

# The output from our simulations is whitespace separated columns. The first 2 columns are the
# simulation ID and the generation. There are then d columns with the population size in each deme, followed
# by d columns with the resistance allele frequency in each deme (where d = number of demes)


## constructor for a simResults class

simResults <- function(x){
  # basic validity checking
  #TODO
  #
  tot.demes = (dim(x)[2] - 2) / 2
  #  Rename the columns
  cnames <- c("id",
              "generation",
              paste("n", 1:tot.demes, sep = "."),
              paste("p", 1:tot.demes, sep = "."))
  names(x) <- cnames
  # rename rows to match generation
  #row.names(x) <- x$generation
  #row.names(x) <- x$generation
  res <- list(num.replicates = length(unique(as.character(x$id))),
              replicate.IDs = unique(as.character(x$id)),
              start.generation = min(x$generation),
              end.generation = max(x$generation),
              sampled.generations = sort(unique(x$generation)),
              num.demes = tot.demes,
              pop.sizes = split(x[,3:(2 + tot.demes)], factor(as.character(x[,1]))),
              allele.freqs = split(x[,(3 + tot.demes): (2 + (tot.demes * 2))], factor(as.character(x[,1]))),
              pop.size.max = max(x[,3:(2 + tot.demes)]), # should be the population size when non-culled
              pop.size.min = min(x[,3:(2 + tot.demes)]) # should be the population size when culled
  )
  #set row names for pop.sizes and allele.freqs to simulation generation
  gens <- as.character(res$sampled.generations)
  #print(gens)
  for (i in 1:length(res$pop.sizes)){
    row.names(res$pop.sizes[[i]]) <- gens[1:length((res$pop.sizes[[i]])[,1])]
  }
  for (i in 1:length(res$allele.freqs))
  {
    row.names(res$allele.freqs[[i]]) <- gens[1:length((res$allele.freqs[[i]])[,1])]
  }
  class(res) <- append(class(res), "simResults")
  return(res)
  
}

## print method for simResults

print.simResults <- function(x){
  cat("simResults object\n")
  cat("number of replicate simulations: ", x$num.replicates, "\n")
  cat("number of demes:", x$num.demes, "\n")
  cat("Duration: from generation", x$start.generation, "to generation", x$end.generation, "\n")
}

## load a simResults object directly from file

read.simResults <- function(file){
  dat <- read.table(file)
  return(simResults(dat))
}

## get the mean allele frequency across demes for a single simulation
## Because simulation where the resistance allele is lost terminate prematurely, by default we pack
## those sims with 0 allele frequencies to fill out the generations that are missing

simulation.mean.allele.freqs <- function(x, simulation.index, pack = T){
  sim <- x$allele.freqs[[simulation.index]]
  mean.freqs <- apply(sim,
                      1,
                      mean)
  if (pack){
    if (length(mean.freqs) < length(x$sampled.generations)){
      mean.freqs <- c(mean.freqs,
                      rep(0.0,
                          times = length(x$sampled.generations) - length(mean.freqs)))
    }
  }
  return(mean.freqs)
}

## get the mean allele frequencies across demes for all simulations
## returns a matrix, rows are sampled generations, columns are replicate sims

mean.allele.freqs <- function(x){
   sim.count = length(x$allele.freqs)
   result <- cbind(simulation.mean.allele.freqs(x, 1))
   for (i in 2: sim.count)
     result <- cbind(result, simulation.mean.allele.freqs(x, i))
   return(result)
                     
}

## get the mean allele frequency over all simulations

grand.mean.allele.freqs <- function(x){
  mean.freqs <- mean.allele.freqs(x)
  grand.mean <- apply(mean.freqs,
                      1,
                      mean)
  return(grand.mean)
}

## get the standard error of the grand mean allele freqs over all simulations

se.grand.mean.allele.freqs <- function(x){
  mean.freqs <- mean.allele.freqs(x)
  
  std.dev <- apply(mean.freqs,
                   1,
                   sd)
  return(std.dev / sqrt(dim(mean.freqs)[2]))
}


## get the mean proportion of demes that are in culled status for a single simulation
## Because simulation where the resistance allele is lost terminate prematurely, by default we pack
## those sims with 0 allele frequencies to fill out the generations that are missing

simulation.cull.proportions <- function(x, simulation.index, pack = T){
  sim <- x$pop.sizes[[simulation.index]]
  culled <- sim == x$pop.size.min
  culled.count <- apply(culled,
                        1,
                        sum)
  if (pack){
    if (length(culled.count) <- length(x$sampled.generations)){
      culled.count <- c(culled.count,
                        rep(0,
                            times = length(x$sampled.generations) - length(culled.count)))
    }
  }
  return(culled.count / x$num.demes)
  
}


#get the proportion of demes that are in culled status for every simulation
# returns a matrix, rows are sampled generations, columns are replicate simulations
cull.proportions <- function(x){
  sim.count <- length(x$pop.sizes)
  result <- cbind(simulation.cull.proportions(x, 1))
  for (i in 2: sim.count){
    result <- cbind(result, simulation.cull.proportions(x, i))
  }
  return(result)
}

# Get the mean proportion of demes in culled status over all simulations

mean.cull.proportions <- function(x){
   cull.props <- cull.proportions(x)
   result <- apply(cull.props,
                    1,
                    mean)
   return(result)
}

# get the standard error of the mean cull proportions

se.mean.cull.proportions <- function(x){
  cull.props <- cull.proportions(x)
  std.dev <- apply(cull.props,
                   1,
                   sd)
  return(std.dev / sqrt(dim(cull.props)[2]))
}


