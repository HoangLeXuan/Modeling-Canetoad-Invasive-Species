

###############################################################
#### Begin of required functions in preparation for main the code

# 1. A dispersal function for calculating what fraction 
# of the population disperses a distance 'd'

dispersal_func <- function(d, Phi){ 
  return(exp(-Phi*d))
}

# 2. Implementation of the Ricker logistic equation 
# This equation dictates the temporal dynamics of each 
# single local population

ricker <- function(xt, r, k, eps, th){
  xt[xt < 0] <- 0
  xt <- xt * exp(r*(1-(xt/k)) + eps)
  xt[xt < th] <- 0
  return(xt)
}

#### End of required functions in preparation for main the code
###############################################################

#### Main application

# this library is needed to perform calculations over the network
# of metapopulation topological structure
require(igraph)

###############################################################
#### Start of model parameters & checking functions 
#### above work fine with the given parameter values

phi <- 1/5              #dispersal distance (5k for rabbits)

##############
# this is how our dispersal kernel looks like...
x <- seq(0,50,0.01)
y <- dispersal_func(x, phi)

plot(x,y,cex=.1, ylab='Fraction of the population moving', xlab='Distance (km)')
###############

r <- 0.4257087          #growth rate
k <- 160                #carrying capacity
sigma <- 0.3149903      #standard deviation of the error
low_th <- 10            #minimum population threshold

################################

################################
#this is how the dynamics of a single population look like:
time_series <- c()
xt <- runif(1, 1 , k)
for(i in 1:100){
  time_series <- append(time_series, xt)
  eps <- rnorm(1, 0, sigma)
  xt <- ricker(xt,r,k,eps,low_th)
}

plot(time_series, cex=.1, type='l')

################################

#### End of model parameters & checking
###############################################################

###############################################################
# Start of setting up the simulation protocol

# Number of replicates to be run for each experiment
replicates <- 20

# This parameter determines the management action to be simulated
# change between 'population' and 'resource' to do population or
# resource reduction respectively
mngmt_action <- 'resource'

dispersal_values <- c(0, 0.3, 0.6, 0.9)          # include all the dispersal values to be tested in this vector
mngmt_extents <- c(.1, .5, .9)                   # include all the management extents (e) to be tested here
mngmt_levels <- c(.3, .6, .9)                    # include all the management levels (l) to be tested here

# management spatial strategies (s) - possible values: random, correlated, hub
mngmt_types <- c('random', 'correlated', 'hub')  

# this is the number of the iteration at which management action is applied
mngmt_iter <- 71

# this is the length (number of iterations) the management is action is applied for
# for the results in the paper: 1 for resource reduction, 9 (which is actually 10 iterations) 
# for population reduction
if(mngmt_action == 'population'){
  mngmt_length <- 9
}else{
  mngmt_length <- 1
}

# number of iterations at which output recording begins and lasts for, respectively
if(mngmt_action == 'population'){
  record_output <- 91
}else{
  record_output <- 73
}
record_length <- 9

# here we calculate the total number of iterations
iterations <- record_output + record_length

# this is where we keep the output
output <- NULL

# vector containing the names of the model metapopulation structures
layouts <- c('star', 'ring', 'neighbours', 'ring-random', 'ring-hub', 'scale-free')

# End of setting up the simulation protocol
###############################################################

###############################################################
# The core of the simulations start here!

# The first loop is thorugh the metapopulation structures (layouts)
# We are going to repeat the same process (replicates of management experiments)
# for each one of the model metapopulation structures

for(grid_id in layouts){
  # here we obtain a network representation of the metapopulation structure
  if(grid_id == 'star'){
    g <- graph.star(15, mode='undirected')
  }else if(grid_id == 'ring'){
    g <- graph.ring(15, directed=FALSE)
  }else if(grid_id == 'neighbours'){
    g <- watts.strogatz.game(1, 15, 2, 0)
  }else if(grid_id == 'scale-free'){
    g <- barabasi.game(15, directed=FALSE)
  }else if(grid_id == 'ring-random'){
    g <- graph.ring(15, directed=FALSE)
    g[from=c(1,3,4,7,9), to=c(10, 12, 15, 11, 5)] <- TRUE
  }else if(grid_id == 'ring-hub'){
    g <- graph.ring(15, directed=FALSE)
    g[1, c(10, 8, 6, 11, 5)] <- TRUE
  }
  
  # we convert the network to a connectivity matrix
  connectivity <- as.matrix(get.adjacency(g))
  
  # and then to a landscape matrix according to the dispersal function
  # This matrix contains information on the fraction of dispersers from 
  # each population going to each other and will be used to simulate dispersal below
  landscape <- connectivity
  landscape <- round(landscape, 1)
  landscape <- dispersal_func(landscape, phi)
  landscape[landscape == 1] <- 0
  landscape <- landscape/rowSums(landscape)
  
  # the number of populations is the number of vertices in the network g
  populations <- vcount(g)
  
  # the carrying capacity of each local population is stored in the carrying vector
  # original_ks is used below to restore k to original values
  original_ks <- carrying <- rep(k,populations)
  
  #### this second loop repeats the management experiments times the number of replicates specified above
  #### using the parameter 'replicates'
  for(repl in 1:replicates){
    print(paste('replicate: ', repl))
    
    ### Initial abundance of each local population (drawn from k +/- 10%)
    init_abs <- rnorm(length(carrying), mean=carrying, sd=(carrying*0.1))
    
    #### the fourth loop goes over the dispersal values (d)
    for(dispersal in dispersal_values){
      # This is the fraction of dispersing individuals
      frac_disp <- rep(dispersal, populations)
      
      #### the fifth loop goes over the management extents (e)
      for(pext in mngmt_extents){
        # This is the number of populations to be managed
        to_perturb <- floor(populations*pext)
        
        #### the sixth loop goes over the spatial management strategies (s)
        for(ptype in mngmt_types){
          # depending on the spatial management strategy chosen, the populations
          # to be managed are different. These are chosen here.
          if(ptype == 'random'){
            perturbed_pops <- sample(populations, to_perturb)   
          }else if(ptype == 'correlated'){
            perturbed_pops <- c()
            chosen <- 0
            while(chosen < to_perturb){
              current <- sample(populations, 1)
              while(current %in% perturbed_pops){
                current <- sample(populations, 1)
              }
              neighs <- igraph::neighbors(g, current)
              perturbed_pops <- append(perturbed_pops, append(current, neighs))
              perturbed_pops <- unique(perturbed_pops)
              chosen <- length(perturbed_pops)
              
            }
            
            if(chosen > to_perturb){
              perturbed_pops <- sample(perturbed_pops, to_perturb)
            }
            
          }else if(ptype == 'hub'){
            hubs <- order(degree(g), decreasing=TRUE)
            perturbed_pops <- hubs[1:to_perturb]
          }
          
          #### finally the last (seventh) loop goes over the spatial management levels (l)
          for(plevel in mngmt_levels){
            ### here we initialise the populations for this replicate
            cur_abs <- init_abs
            dynamics <- as.matrix(t(cur_abs))
            carrying <- original_ks
            
            # this is where the actual simulation of the population dynamics occur
            for(i in 1:iterations){         
              # here we check if we are within the period of management
              # if we are, then apply management
              if(i >= mngmt_iter & i < (mngmt_iter + mngmt_length)){
                
                if(mngmt_action == 'population'){
                  #this is for poison baiting (population reduction)
                  cur_abs[perturbed_pops] <- cur_abs[perturbed_pops] - (cur_abs[perturbed_pops]*plevel)
                }else if(mngmt_action == 'resource'){
                  #this is for warren ripping (affecting k)
                  carrying[perturbed_pops] <- original_ks[perturbed_pops] - (original_ks[perturbed_pops]*plevel)
                  above <- which(cur_abs[perturbed_pops] > carrying[perturbed_pops])
                  cur_abs[above] <- carrying[above]
                }
              }
              
              # stochasticity for population growth
              eps <- rnorm(populations, 0, sigma)
              
              # here populations grow according to ricker logistic function
              cur_abs <- ricker(cur_abs, r, carrying, eps, low_th)
              dynamics <- rbind(dynamics, cur_abs)        
              
              # here dispersal is simulated by sending a fraction of each local populations
              # to the respective connected populations
              dispersing <- cur_abs * frac_disp
              cur_abs <- cur_abs - dispersing
              arriving <- dispersing*landscape
              arriving <- colSums(arriving)
              
              cur_abs <- cur_abs + arriving
              
            }
            
            #### After performing the execution of the dynamics we obtain summary statistics for the
            #### resulting dynamics
            # these populations became extinct
            ext_pops <- which(apply(dynamics[(mngmt_iter-10):mngmt_iter,], 2, mean) == 0)
            
            if(length(ext_pops) > 0){
              dynamics <- dynamics[,-ext_pops]
            }
            
            if(is.null(dim(dynamics))){
              dyns_b4 <- dynamics[1:mngmt_iter]
            }else{
              dyns_b4 <- dynamics[1:mngmt_iter,]
            }
            
            # here we calculate the coefficient of variation
            if(is.null(dim(dyns_b4))){
              if(length(dyns_b4) > 0){
                mean_cv <- cv_reg <- sd(dyns_b4)/mean(dyns_b4)
              }else{
                mean_cv <- 0
                cv_reg <- 0
              }
            }else{
              mean_cv <- mean(apply(dyns_b4, 2, sd)/apply(dyns_b4, 2, mean))
              regional_dyns <- rowSums(dyns_b4)
              cv_reg <- sd(regional_dyns)/mean(regional_dyns)
            }
            
            if(is.null(dim(dynamics))){
              ext_pops_during <- which(mean(dynamics[mngmt_iter:(mngmt_iter+mngmt_length)]) == 0)
              
            }else{
              ext_pops_during <- which(apply(dynamics[mngmt_iter:(mngmt_iter+mngmt_length),], 2, mean) == 0)
              if(length(ext_pops_during) > 0){
                dynamics <- dynamics[,-ext_pops_during]
              }
            }
            
            
            # here we calculate the mean, max and min abundaces during management for the
            # period specified by the parameters above
            if(is.null(dim(dynamics))){
              window <- dynamics[(mngmt_iter-10):(mngmt_iter)]
            }else{
              window <- dynamics[(mngmt_iter-10):(mngmt_iter),]
            }
            
            if(is.null(dim(window))){
              if(length(window) > 0){
                mean_10_during <- mean(window)
                max_10_during <- max(window)
                min_10_during <- min(window) 
              }else{
                mean_10_during <- 0
                max_10_during <- 0
                min_10_during <- 0
              }
            }else{
              regional <- rowSums(window)
              mean_10_during <- mean(regional)
              max_10_during <- max(regional)
              min_10_during <- min(regional)
            }
            
            # here we calculate the mean, max and min abundaces after management for the
            # period specified by the parameters above
            if(is.null(dim(dynamics))){
              window <- dynamics[record_output:(record_output+record_length)]
              ext_pops_after <- which(mean(window) == 0)
              if(length(ext_pops_after) > 0){
                window <- c()
              }
            }else{
              window <- dynamics[record_output:(record_output+record_length),]
              ext_pops_after <- which(apply(window, 2, mean) == 0)
              if(length(ext_pops_after) > 0){
                window <- window[,-ext_pops_after]
              }
            }
            
            if(is.null(dim(window))){
              if(length(window) > 0){
                
                mean_10_after <- mean(window)
                max_10_after <- max(window)
                min_10_after <- min(window)
                
                mean_cv_10_after <- cv_reg_10_after <- sd(window)/mean(window)
              }else{
                mean_10_after <- 0
                max_10_after <- 0
                min_10_after <- 0
                
                mean_cv_10_after <- 0
                cv_reg_10_after <- 0
              }
            }else{
              regional <- rowSums(window)
              mean_10_after <- mean(regional)
              max_10_after <- max(regional)
              min_10_after <- min(regional)
              
              mean_cv_10_after <- mean(apply(window, 2, sd)/apply(window, 2, mean))
              cv_reg_10_after <- sd(regional)/mean(regional)
            }
            
            ext_pops <- length(ext_pops)
            ext_pops_after <- length(ext_pops_after)
            cur_res <- data.frame(grid_id, repl, dispersal,
                                  plevel=as.character(plevel), ptype, pext=as.character(pext),
                                  ext_pops, ext_pops_after, mean_cv, cv_reg,
                                  mean_10_during, min_10_during, max_10_during,
                                  mean_10_after, max_10_after, min_10_after,
                                  mean_cv_10_after, cv_reg_10_after)
            
            if(is.null(output)){
              output <- cur_res
            }else{
              output <- rbind(output, cur_res)
            }
          }
        }
      }
    }
  }
}

# End of simulation experiments
###############################################################


##########################################################################
# Start of saving/reporting results

# Executing this instruction saves the output to a file in csv format
write.csv(output, file=paste('./model-landscapes-results-',mngmt_action,'.csv',sep=''))

############
# Plotting
# We can also plot the results
# Execute the code below to reproduce the panels in Figure 1 in the manuscript
# Pdf's will be saved to your working directory

# Note: We use the GGPlot library for plotting

require(grid)
require(gridBase)
require(ggplot2)

# once you have ran the experiments for both management actions
# assign the correspoding outputs to the variables below
# so the plots will contain both results

# output for population reduction
output_pops <- output

# output for resource reduction
output_res <- output

plots <- list()
idx <- 1

plot_features <- "
stat_smooth(method='loess', se=TRUE, size = 1, na.rm=TRUE) +
theme(panel.background=element_blank(), plot.margin=unit(c(1,.5,-.5,0),'cm'),
legend.position=c(-.1, 1.15), legend.direction='horizontal', 
axis.title.y = element_text(size=20), axis.title.x = element_text(size=20)) + "

remove_legend <- "scale_colour_discrete(guide=FALSE) "

for(l in layouts){
  
  if(l == 'star'){
    g <- graph.star(15, mode='undirected')
  }else if(l == 'ring'){
    g <- graph.ring(15, directed=FALSE)
  }else if(l == 'neighbours'){
    g <- watts.strogatz.game(1, 15, 2, 0)
  }else if(l == 'scale-free'){
    g <- barabasi.game(15, directed=FALSE)
  }else if(l == 'ring-random'){
    g <- graph.ring(15, directed=FALSE)
    g[from=c(1,3,4,7,9), to=c(10, 12, 15, 11, 5)] <- TRUE
  }else if(l == 'ring-hub'){
    g <- graph.ring(15, directed=FALSE)
    g[1, c(10, 8, 6, 11, 5)] <- TRUE
  }
  
  data <- output_pops
  data$dispersal <- as.numeric(data$dispersal)
  data$ext_pops_after <- as.numeric(data$ext_pops_after)
  data$plevel <- as.numeric(as.character(data$plevel))
  
  temp <- data[which(data$dispersal == 0.3  & data$plevel == .6 & data$grid_id == l),]
  temp$ptype <- as.factor(temp$ptype)
  temp$ptype <- relevel(temp$ptype, 'random')
  
  temp$pext <- as.numeric(as.character(temp$pext))
  
  if(l == 'star'){
    ylab <- " + ylab('Max Abundance')"
  }else{
    ylab <- " + ylab('')"
  }
  
  xlab <- " + xlab('')"
  
  axis_limits <- " coord_cartesian(ylim = c(.6,1.05), xlim=c(0,1)) + "
  
  p2 <- eval(parse(text=paste(
    "ggplot(temp, aes(pext, ((max_10_after/max_10_during)), colour=as.factor(ptype))) +",
    plot_features, axis_limits, remove_legend, ylab, xlab)))
  
  data <- output_res
  data$dispersal <- as.numeric(data$dispersal)
  data$ext_pops_after <- as.numeric(data$ext_pops_after)
  data$plevel <- as.numeric(as.character(data$plevel))
  
  temp <- data[which(data$dispersal == 0.3  & data$plevel == .6 & data$grid_id == l),]
  temp$ptype <- as.factor(temp$ptype)
  temp$ptype <- relevel(temp$ptype, 'random')
  
  temp$pext <- as.numeric(as.character(temp$pext))
  
  xlab <- " + xlab('Management Extent')"
  axis_limits <- " coord_cartesian(ylim = c(.4,1.05), xlim=c(0,1)) + "
  
  p3 <- eval(parse(text=paste(
    "ggplot(temp, aes(pext, ((max_10_after/max_10_during)), colour=as.factor(ptype))) +",
    plot_features, axis_limits, remove_legend, ylab, xlab)))
  
  pdf(paste(l, '-test.pdf',sep=''),width=3,height=8)
  plot.new()
  vp1 <- viewport(x=0, y=0.66, width=1, height=0.33, just = c("left", "bottom"))
  vp2 <- viewport(x=0, y=0.21, width=1, height=0.6, just = c("left", "bottom"))
  vp3 <- viewport(x=0, y=0.02, width=1, height=0.6, just = c("left", "bottom"))
  
  pushViewport(vp3)
  print(p3, vp = vp3)
  
  upViewport()
  pushViewport(vp2)
  print(p2, vp = vp2)
  
  upViewport()
  pushViewport(vp1)
  par(new=TRUE, fig=gridFIG(), mar=c(0,0,0,0))
  
  if(l=='star' | l=='scale-free' | l=='neighbours'){
    plot(g, vertex.color=NA, vertex.size=22)
  }else{
    plot(g, vertex.color=NA, layout=layout.circle, vertex.size=22)  
  }
  
  dev.off()
  
}


