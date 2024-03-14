


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

# 3. This function calculates the geodesic distance between 
# two points specified by radian latitude/longitude using the
# Spherical Law of Cosines (slc)
# This function was obtained here:
# http://www.r-bloggers.com/great-circle-distance-calculations-in-r/
gcd.slc <- function(long1, lat1, long2, lat2) {
  R <- 6371 # Earth mean radius [km]
  d <- acos(sin(lat1)*sin(lat2) + cos(lat1)*cos(lat2) * cos(long2-long1)) * R
  
  if(is.nan(d)){
    print(long1 != long2)
    print(lat1 != lat2)
    print(paste(long1, lat1, long2, lat2))
  }
  
  return(d) # Distance in km
}

# 4. This function converts degrees to radians
deg2rad <- function(deg) return(deg*pi/180)

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

# this is the threshold distance (in km) at which local populations 
# are considered to be connected. i.e., any to populations this distance
# or less away are connected on the metapopulation structure
link_distance <- 10     

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

# End of setting up the simulation protocol
###############################################################

###############################################################
# The core of the simulations start here!

# The first loop is thorugh the real-world landscapes
# We are going to repeat the same process (replicates of management experiments)
# for each one of the real-world metapopulations

# there are 10 such landscapes available
for(grid_id in 1:10){
  print('grid id')
  print(grid_id)
  
  # change the following line if you want to read in other landscapes.
  # the filename of the landscape should be specified here
  data <- read.csv(paste('./landscapes/tbl_MetapopGrid1km_landsc50km_', grid_id, '.csv',sep=''), header=T, sep=',')
  pops <- which(data$PopSize > 0)
  
  if(length(pops) == 0){
    next
  }
  
  # the carrying capacity of each local population is stored in the carrying vector
  # original_ks is used below to restore k to original values
  original_ks <- carrying <- (data$PopSize[pops])*100
  
  # we obtain the connectivity matrix calculating the distances between
  # local cells in the real-world landscapes
  connectivity <- matrix(0, nrow=length(pops), ncol=length(pops))
  
  idx_j <- 0
  for(j in pops){
    idx_j <- idx_j + 1
    x_j <- deg2rad(data[j,]$XCoordCell)
    y_j <- deg2rad(data[j,]$YCoordCell)
    
    idx_k <- 0
    for(pop in pops){
      idx_k <- idx_k + 1
      x_k <- deg2rad(data[pop,]$XCoordCell)
      y_k <- deg2rad(data[pop,]$YCoordCell)
      
      if(x_j != x_k & y_j != y_k){
        d_jk <- gcd.slc(x_j, y_j, x_k, y_k)
        if(d_jk <= link_distance){
          connectivity[idx_j,idx_k] <- d_jk
        }
      }
    }
  }
  
  g <- graph.adjacency(connectivity, mode='undirected', weighted=TRUE)
  
  #we remove the isolated populations from the landscape matrix and from the graph
  isolates <- which(igraph::degree(g) == 0)
  if(length(isolates > 0)){
    g <- delete.vertices(g, isolates)
    connectivity <- connectivity[-isolates, -isolates]
  }
  
  # and then to a landscape matrix according to the dispersal function
  # This matrix contains information on the fraction of dispersers from 
  # each population going to each other and will be used to simulate dispersal below
  landscape <- connectivity
  landscape <- round(landscape, 1)
  landscape <- dispersal_func(landscape, phi)
  landscape[landscape == 1] <- 0
  landscape <- landscape/rowSums(landscape)
  
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
                                  ext_pops, ext_pops_after, frac_living_pops=((pops - ext_pops_after)/pops),
                                  mean_cv, cv_reg,
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
write.csv(output, file=paste('./real-world-landscapes-results-',mngmt_action,'.csv',sep=''))

############
# Plotting
# We can also plot the results
# Execute the code below to reproduce Figure 3 in the manuscript
# Pdf's will be saved to your working directory

# Note: We use the GGPlot library for plotting

require(ggplot2)
require(grid)

vplayout <- function(x, y) viewport(layout.pos.row=x, layout.pos.col=y)

temp <- output

temp$pext <- as.numeric(as.character(temp$pext))

plot_features <- "
facet_wrap(~dispersal, scales='fixed', nrow=1) + 
stat_smooth(method='loess', method.args=list(degree=1), se=TRUE, size = 1, na.rm=TRUE) +
theme(panel.background=element_blank(), plot.margin=unit(c(0,1,0,0),'cm'),
legend.position=c(.5,.82), legend.direction='horizontal', 
legend.title=element_text(size=20), legend.text=element_text(size=15),
axis.title.y = element_text(size=20), axis.title.x = element_text(size=20),
plot.title = element_text(size = 20), strip.background = element_blank(), 
strip.text = element_text(size = 20)) + 
xlab(' Management Extent       Management Extent       Management Extent       Management Extent ') +"

remove_legend <- "scale_colour_discrete(guide=FALSE) + "

p1 <- eval(parse(text=paste(
  "ggplot(temp, aes(x=pext, y=((cv_reg_10_after*100)/(cv_reg*100)), color=as.factor(plevel))) +",
  plot_features,
  "scale_colour_discrete(name='Management level') +  ylab('CV') + xlab('') +
  ggtitle('Dispersal')")))

p2 <- eval(parse(text=paste(
  "ggplot(temp, aes(pext, max_10_after/max_10_during, color=as.factor(plevel))) +",
  plot_features, remove_legend,
  "ylab('Max Abundance') + xlab('') + theme(strip.background = element_blank(), strip.text = element_blank())")))

p3 <- eval(parse(text=paste(
  "ggplot(temp, aes(pext, frac_living_pops, color=as.factor(plevel))) +",
  plot_features, remove_legend,
  "ylab('Surviving Populations') + 
  theme(strip.background = element_blank(), strip.text = element_blank(),
  plot.margin=unit(c(0,1,0.1,0),'cm'))")))


# Change here the name of the output file
pdf(paste('test.pdf',sep=''),width=13,height=9)

grid.newpage()
pushViewport(viewport(layout=grid.layout(3,1)))

print(p1,vp=vplayout(1,1))
print(p2,vp=vplayout(2,1))
print(p3,vp=vplayout(3,1))

dev.off()
