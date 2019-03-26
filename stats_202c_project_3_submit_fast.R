require(grDevices)
library(ggplot2)
library(knitr)
library(rstackdeque)

options("expressions"=20000) 
options("max-ppsize"=500000)

findNeighbors = function(x,y,nrow,ncol) {
    
    # can only have a max of 4 neighbors 
    neighbors = matrix(0, nrow=4, ncol=2)
    
    # check top
    if ((y-1) > 0) {
        neighbors[1,2] = x
        neighbors[1,1] = y-1      
    } else { # wrap around
        neighbors[1,2] = x
        neighbors[1,1] = nrow
    }
    
    # check bottom
    if ((y+1) <= nrow) {
        neighbors[2,2] = x
        neighbors[2,1] = y+1
    } else { # wrap around
        neighbors[2,2] = x
        neighbors[2,1] = 1
    }
    
    # check left
    if ((x-1) > 0) {
        neighbors[3,2] = x-1
        neighbors[3,1] = y
    } else { # wrap around
        neighbors[3,2] = ncol
        neighbors[3,1] = y
    }
    
    # check right
    if ((x+1) <= ncol) {
        neighbors[4,2] = x+1
        neighbors[4,1] = y
    } else { # wrap around
        neighbors[4,2] = 1
        neighbors[4,1] = y
    }
    
    return(neighbors)
}

turnEdgesOnOff = function(mc, qe) {
    n=nrow(mc)
    iedges_mc = matrix(0,n,n)
    jedges_mc = matrix(0,n,n)
    
    for (ii in 1:n) {
        for (jj in 1:n) {
            
            # edge in the i direction
            if (ii == n) {
                inext = 1 # wrap around
            } else {
                inext = ii + 1
            }
            if ((mc[ii,jj] == mc[inext,jj]) && (runif(1) < qe)) {
                iedges_mc[ii,jj] = 1
            }
            
            # edge in the j direction
            if (jj == n) {
                jnext = 1 # wrap around
            } else {
                jnext = jj + 1
            }
            if ((mc[ii,jj] == mc[ii,jnext]) && (runif(1) < qe)) {
                jedges_mc[ii,jj] = 1
            }
        }
    }
    
    return(list(iedges = iedges_mc, jedges = jedges_mc))
}

properLabel = function(labelLabel, label) {
    label_ = label
    while(labelLabel[label_] != label_) {
        label_ = labelLabel[label_]
    }
    return(label_)
}

labelClusters = function(iedges_mc, jedges_mc){
    
    labelLabel = matrix(0,nrow(iedges_mc)*ncol(iedges_mc),1)
    cluster = matrix(0, nrow(iedges_mc), ncol(iedges_mc))
    label = 1
    for (ii in 1:n) {
        for (jj in 1:n) {
            
            # find previously visited sites connected to i, j
            bonds = 0
            iBond = matrix(0,4,1)
            jBond = matrix(0,4,1)
            
            # check bond to i-1,j
            if (ii > 1 && iedges_mc[ii-1,jj]) {
                bonds = bonds+1
                iBond[bonds] = ii-1
                jBond[bonds] = jj
            }
            
            # wrap around at the boundary (if i,j is the last site, check
            # bond to i+1,j)
            if (ii == n && iedges_mc[ii,jj]) {
                bonds = bonds+1
                iBond[bonds] = 1
                jBond[bonds] = jj
            }
            
            # check bond to i,j-1
            if (jj > 1 && jedges_mc[ii,jj-1]) {
                bonds = bonds+1
                iBond[bonds] = ii
                jBond[bonds] = jj-1
            }
            
            # wrap around at the bondary at the last site
            if (jj == n && jedges_mc[ii,jj]) {
                bonds = bonds+1
                iBond[bonds] = ii
                jBond[bonds] = 1
            }
            
            # check number of bonds to previously visited sites
            if (bonds == 0) { 
                # need to start a new cluster
                cluster[ii,jj] = label
                labelLabel[label] = label
                label = label + 1                
            } else {
                # re-label bonded sites with smallest proper label
                minLabel = label
                for (b in 1:bonds) {
                    pLabel = properLabel(labelLabel, cluster[iBond[b],jBond[b]])
                    if (minLabel > pLabel) {
                        minLabel = pLabel
                    }
                }          
                
                # set current site label to smallest proper label
                cluster[ii,jj] = minLabel
                
                # reset the proper label links on the previous labels
                for (b in 1:bonds) {
                    pLabel = cluster[iBond[b],jBond[b]]
                    labelLabel[pLabel] = minLabel
                    
                    # reset label on connected sites
                    cluster[iBond[b],jBond[b]] = minLabel
                }
            }
        }
    }
    
    return(list(cluster=cluster, labelLabel=labelLabel))
}


flipClusterSpins_V1 = function(cluster_in, mc, labelLabel) {
    n = nrow(cluster_in)
    cluster = cluster_in
    new_mc = mc
    flips = 0 # count number of spins that are flipped
    
    # keep track of new cluster flipped
    sNewChosen = matrix(0,n*n,1)
    # new cluster spin values
    sNew = matrix(0,n*n,1)
    
    for (ii in 1:n) {
        for (jj in 1:n) {
            # random new cluster spins values have not been set
            nn = (ii-0)*n+jj
            sNewChosen[nn] = 0
            
            # replace all labels by their proper values
            cluster[ii,jj] = properLabel(labelLabel, cluster[ii,jj])           
        }
    }
    
    for (ii in 1:n) {
        for (jj in 1:n) {
            # find the new proper label of the cluster
            label = cluster[ii,jj]
            
            # choose a random new spin value for cluster only
            # if this has not already been done
            if (sNewChosen[label] == 0) {
                rand = runif(1)
                if (rand < 0.5 ) {
                    sNew[label] = 1
                } else {
                    sNew[label] = 0
                }   
                sNewChosen[label] = 1
            }
            
            # reset the spin value and count number of flips
            if (new_mc[ii,jj] != sNew[label]) {
                new_mc[ii,jj] = sNew[label]
                flips = flips + 1
            }
        }
    }
    
    # get average size of CP for this iteration
    cluster_labels = unique(matrix(cluster, n*n,1))
    CP_size = matrix(0, length(cluster_labels),1)
    ll=1
    for (cl in cluster_labels) {
        CP_size[ll] = length(which(cluster==cl))
        ll=ll+1
    }
    
    return(list(cluster=cluster, new_mc=new_mc, ave_CP_size=mean(CP_size)))
}

flipClusterSpins_V2 = function(cluster_in, mc, labelLabel, rr, cc) {
    n = nrow(cluster_in)
    cluster = cluster_in
    new_mc = mc
    flips = 0 # count number of spins that are flipped
    
    # keep track of new cluster flipped
    sNewChosen = matrix(0,n*n,1)
    # new cluster spin values
    sNew = matrix(0,n*n,1)
    
    for (ii in 1:n) {
        for (jj in 1:n) {
            # random new cluster spins values have not been set
            nn = (ii-0)*n+jj
            sNewChosen[nn] = 0
            
            # replace all labels by their proper values
            cluster[ii,jj] = properLabel(labelLabel, cluster[ii,jj])           
        }
    }
    
    # randomly choose a seed from the grid
    rr=round(runif(1, min=1, max=n))
    cc=round(runif(1, min=1, max=n))
    
    # flip the state 
    new_state = as.numeric(!mc[rr,cc])
    
    # find the new proper label of the cluster
    label = cluster[rr,cc]
    
    # flip the entire cluster
    new_mc[cluster==label] = new_state
    
    return(list(cluster=cluster, new_mc=new_mc, ave_CP_size=sum(cluster==label)))
}

up = function(rr, cc, n) {
    coords=matrix(0,2,1)
    if (rr == 1) new_rr = n else new_rr=rr-1
    coords[1]=new_rr
    coords[2]=cc
    return(coords)
}
down = function(rr, cc, n) {
    coords=matrix(0,2,1)
    if (rr == n) new_rr = 1 else new_rr=rr+1
    coords[1]=new_rr
    coords[2]=cc
    return(coords)
}
left = function(rr, cc, n) {
    coords=matrix(0,2,1)
    if (cc == 1) new_cc = n else new_cc=cc-1
    coords[1]=rr
    coords[2]=new_cc
    return(coords)
}
right = function(rr, cc, n) {
    coords=matrix(0,2,1)
    if (cc == n) new_cc = 1 else new_cc=cc+1
    coords[1]=rr
    coords[2]=new_cc
    return(coords)
}

# Wolff algorithm 
growCluster = function(i, j, clusterSpin, mc, cluster, qe) {
    new_cluster = cluster
    new_mc = mc
    n=nrow(mc)
    
    # mark as belonging to the cluster
    new_cluster[i,j] = 1
    # flip it
    new_mc[i,j] = as.numeric(!mc[i,j])
    
    # find the indices of the 4 neighbors
    # assuming periodic boundary conditions
    if(i == 1) iPrev = n else iPrev = i-1
    if(i == n) iNext = 1 else iNext = i+1
    if(j == 1) jPrev = n else jPrev = j-1
    if(j == n) jNext = 1 else jNext = j+1
    
    # if the neighbor spin does not belong to the
    # cluster, then try to add it to the cluster
    if (!new_cluster[iPrev, j]) {
        temp = tryAdd(iPrev, j, clusterSpin, new_mc, qe, new_cluster);
        new_mc = temp$new_mc
        new_cluster = temp$new_cluster
    }
    if (!new_cluster[iNext, j]) {
        temp = tryAdd(iNext, j, clusterSpin, new_mc, qe, new_cluster);
        new_mc = temp$new_mc
        new_cluster = temp$new_cluster
    }
    if (!new_cluster[i, jPrev]) {
        temp = tryAdd(i, jPrev, clusterSpin, new_mc, qe, new_cluster);
        new_mc = temp$new_mc
        new_cluster = temp$new_cluster
    }
    if (!new_cluster[i, jNext]) {
        temp = tryAdd(i, jNext, clusterSpin, new_mc, qe, new_cluster); 
        new_mc = temp$new_mc
        new_cluster = temp$new_cluster
    }
    
    #     image(new_mc*256)
    
    return(list(new_mc=new_mc, new_cluster=new_cluster))
}

tryAdd = function(i, j, clusterSpin, mc, qe, cluster) {
    new_mc = mc
    new_cluster = cluster
    
    if (mc[i,j] == clusterSpin) {
        if (runif(1) < qe) {
            temp = growCluster(i, j, clusterSpin, mc, cluster, qe)
            new_mc = temp$new_mc
            new_cluster = temp$new_cluster
        }
    }
    return(list(new_mc=new_mc, new_cluster=new_cluster))
}

# Wolff algorithm (wrong...stack is too big)
growOneCpAndFlip = function(mc, qe, rr, cc) {
    n=nrow(mc) 
    s = rstack() # store the path so we can backtrack
    cp = matrix(0, nrow=nrow(mc), ncol=ncol(mc)) # store the sites for this CP
    visited = matrix(0, nrow=nrow(mc), ncol=ncol(mc))
    new_mc = mc
    
    # decide the new spin for this CP that we are growing
    # if (runif(1) < 0.5) new_state = 0 else new_state = 1
    
    #     # randomly pick a pixel to grow a CP from
    #     rr=round(runif(1, min=1, max=n))
    #     cc=round(runif(1, min=1, max=n))
    
    # flip the state 
    new_state = as.numeric(!mc[rr,cc])
    
    # save the original state
    original_state = mc[rr,cc]
    
    # insert the current position
    s = insert_top(s, c(rr,cc))
    cp[rr,cc] = 1
    new_mc[rr,cc] = new_state
    visited[rr,cc] = 1
    
    while(!empty(s)) {
        
        coords=peek_top(s)
        s = without_top(s)        
        
        neighbors = findNeighbors(coords[1],coords[2], n, n)
        neighbors = neighbors[cp[neighbors] == 0,]                              
        #         if(sum(cp[neighbors]) == 4) {
        #             next
        #         }
        neighbors = matrix(neighbors, ncol=2)
        if (nrow(neighbors) == 0) {
            next
        }
        
        same_state = neighbors[mc[neighbors] == original_state,]
        same_state = matrix(same_state, ncol=2)
        if (nrow(same_state) == 0) {
            next
        }
        
        add_to_cluster = same_state[runif(nrow(same_state)) < qe ,]
        add_to_cluster = matrix(add_to_cluster, ncol=2)
        if (nrow(add_to_cluster) == 0) {
            next
        }
        cp[add_to_cluster] = 1
        visited[add_to_cluster] = 1
        new_mc[add_to_cluster] = new_state
        for (aa in 1:nrow(add_to_cluster)) {
            s = insert_top(s, add_to_cluster[aa,])
        } 
        
        if (all(visited)) {
            break
        }
        
        #         temp = mc[neighbors] == original_state
        #         same_state = neighbors[temp,]
        #         if (sum(temp) == 0) {
        #             next
        #         } else if (sum(temp) == 1) {
        #             if (runif(1) < qe) {
        #                 cp[same_state[1], same_state[2]] = 1
        #                 new_mc[same_state[1], same_state[2]] = new_state
        #                 s = insert_top(s, same_state)
        #             }
        #         } else {
        #             
        #             add_to_cluster = same_state[runif(nrow(same_state)) < qe ,]
        #             cp[add_to_cluster] = 1
        #             new_mc[add_to_cluster] = new_state
        #             
        #             if (length(add_to_cluster)==0) {
        #                 next
        #             } else if (length(add_to_cluster)==2) { # only 1 coord
        #                 s = insert_top(s, add_to_cluster)
        #             } else {               
        #                 for (aa in 1:nrow(add_to_cluster)) {
        #                     s = insert_top(s, add_to_cluster[aa,])
        #                 } 
        #             }
        #         }
        
        
        
        #         # traverse order is up, right, down, left (with wrap around)
        #         if (runif(1) < qe) {
        #             coords = up(coords[1],coords[2],n)
        #             if(mc[coords[1],coords[2]] == original_state) {
        #                 s = insert_top(s, coords)
        #                 cp[coords[1],coords[2]] = 1
        #                 new_mc[coords[1], coords[2]] = new_state
        #             }
        #         }     
        #         if (runif(1) < qe) {
        #             coords = right(coords[1],coords[2],n)
        #             if(mc[coords[1],coords[2]] == original_state) {
        #                 s = insert_top(s, coords)
        #                 cp[coords[1],coords[2]] = 1
        #                 new_mc[coords[1], coords[2]] = new_state
        #             }
        #         }
        #         if (runif(1) < qe) {
        #             coords = down(coords[1],coords[2],n)
        #             if(mc[coords[1],coords[2]] == original_state) {
        #                 s = insert_top(s, coords)
        #                 cp[coords[1],coords[2]] = 1
        #                 new_mc[coords[1], coords[2]] = new_state
        #             }
        #         }
        #         if (runif(1) < qe) {
        #             coords = left(coords[1],coords[2],n)
        #             if(mc[coords[1],coords[2]] == original_state) {
        #                 s = insert_top(s, coords)
        #                 cp[coords[1],coords[2]] = 1
        #                 new_mc[coords[1], coords[2]] = new_state
        #             }
        #         }        
    }
    
    #         # debug
    #         image(new_mc*256)
    
    return(list(new_mc=new_mc, CP_size=sum(cp)))
}


#############
# Test case
#############
set.seed(1234567)
n=64
nrow = n
ncol = n
maxIter = 350
betas = c(0.65, 0.75, 0.85, 1.0)
diff_eps = 1e-3


for (sw in 1:2) {
    
    taus = matrix(0, nrow=length(betas), ncol=1)
    h_mc1 = matrix(0,maxIter,1)
    h_mc2 = matrix(0,maxIter,1)
    h_at_convergence = matrix(0,nrow=length(betas), ncol=1)
    beta_CP_size_mc1 = matrix(0,nrow=length(betas), ncol=1)
    beta_CP_size_mc2 = matrix(0,nrow=length(betas), ncol=1)
    bb = 1
    
    if (sw==1) algo="Swendsen-Wang" else algo="Wolff"
    
    for (beta in betas) {
        
        # print(paste0("beta = ", beta))
        
        ave_CP_size_mc1 = matrix(0,nrow=maxIter, ncol=1)
        ave_CP_size_mc2 = matrix(0,nrow=maxIter, ncol=1)
        
        # re-initialize mc1 and mc2
        mc1 = matrix(1, nrow=n, ncol=n)
        # checkerboard
        temp = matrix(0,nrow=n,ncol=1)
        temp[seq(1,n,2)] = 1
        mc2=matrix(0,nrow=n,ncol=n)
        mc2[,seq(1,n,2)] = temp
        temp = matrix(0,nrow=n,ncol=1)
        temp[seq(2,n,2)] = 1
        mc2[,seq(2,n,2)] = temp
        
        for (it in 1:maxIter) {
            
            # print(paste0("iteration: ", it))
            
            # probability of forming an edge, given xs == xt
            qe = 1-exp(-beta)
            
            if (sw == 1) { 
                
                # Swendsen-Wang (form CP's over the entire image and
                # randomly set the spins of each cluster to 0 or 1)
                
                # update the edges for this sweep
                edges = turnEdgesOnOff(mc1, qe)
                iedges_mc1 = edges$iedges
                jedges_mc1 = edges$jedges
                edges = turnEdgesOnOff(mc2, qe)
                iedges_mc2 = edges$iedges
                jedges_mc2 = edges$jedges
                
                # use Hoshen-Kopelman algorithm to identify and label clusters/cps
                temp = labelClusters(iedges_mc1, jedges_mc1)
                cluster_mc1 = temp$cluster
                labelLabel_mc1 = temp$labelLabe
                temp = labelClusters(iedges_mc2, jedges_mc2)
                cluster_mc2 = temp$cluster
                labelLabel_mc2 = temp$labelLabel
                
                # flip the clusters randomly (ber(0.5))
                temp = flipClusterSpins_V1(cluster_mc1, mc1, labelLabel_mc1)
                cluster_mc1 = temp$cluster
                new_mc1 = temp$new_mc
                ave_CP_size_mc1[it] = temp$ave_CP_size
                temp = flipClusterSpins_V1(cluster_mc2, mc2, labelLabel_mc2)
                cluster_mc2 = temp$cluster
                new_mc2 = temp$new_mc
                ave_CP_size_mc2[it] = temp$ave_CP_size
                
            } else {
                
                if (beta < 1.0) {
                    # Version 2 (Wolff algorithm) (traversing the grid for a single 
                    # sweep instead of randomly picking seeds)
                    # uses recursion (seems to be faster, but can't handle beta = 1.0)
                    for (rr in 1:n) {
                        for (cc in 1:n) {
                            cluster_mc1 = matrix(0, n, n)
                            cluster_mc2 = matrix(0, n, n)
                            
                            temp = growCluster(rr, cc, mc1[rr, cc], mc1, cluster_mc1, qe) 
                            new_mc1 = temp$new_mc
                            cluster_mc1 = temp$new_cluster
                            temp = growCluster(rr, cc, mc2[rr, cc], mc2, cluster_mc2, qe)
                            new_mc2 = temp$new_mc
                            cluster_mc2 = temp$new_cluster
                            
                            #                         # debug
                            #                         image(new_mc1*256)
                            #                         image(new_mc2*256)
                            
                            ave_CP_size_mc1[it] = ave_CP_size_mc1[it] + sum(cluster_mc1)
                            ave_CP_size_mc2[it] = ave_CP_size_mc2[it] + sum(cluster_mc2)
                            
                            mc1 = new_mc1
                            mc2 = new_mc2
                        }
                    }
                } else {
                    
                    #                     # wolff
                    #                     # uses a stack
                    #                     for (rr in 1:n) {
                    #                         for (cc in 1:n) {
                    #                             temp = growOneCpAndFlip(mc1, qe, rr, cc)
                    #                             new_mc1 = temp$new_mc
                    #                             ave_CP_size_mc1[it] = ave_CP_size_mc1[it] + temp$CP_size
                    #                             temp = growOneCpAndFlip(mc2, qe, rr, cc)
                    #                             new_mc2 = temp$new_mc
                    #                             ave_CP_size_mc2[it] = ave_CP_size_mc2[it] + temp$CP_size
                    #                             
                    #                             mc1 = new_mc1
                    #                             mc2 = new_mc2
                    #                         }
                    #                     }
                    
                    
                    # Wolff
                    # do it like how we did SW1, but pick a random cluster and flip the sites
                    # update the edges for this sweep
                    for (rr in 1:n) {
                        for (cc in 1:n) {
                            edges = turnEdgesOnOff(mc1, qe)
                            iedges_mc1 = edges$iedges
                            jedges_mc1 = edges$jedges
                            edges = turnEdgesOnOff(mc2, qe)
                            iedges_mc2 = edges$iedges
                            jedges_mc2 = edges$jedges
                            
                            # use Hoshen-Kopelman algorithm to identify and label clusters/cps
                            temp = labelClusters(iedges_mc1, jedges_mc1)
                            cluster_mc1 = temp$cluster
                            labelLabel_mc1 = temp$labelLabe
                            temp = labelClusters(iedges_mc2, jedges_mc2)
                            cluster_mc2 = temp$cluster
                            labelLabel_mc2 = temp$labelLabel
                            
                            # flip the chosen cluster
                            temp = flipClusterSpins_V2(cluster_mc1, mc1, labelLabel_mc1, rr, cc)
                            cluster_mc1 = temp$cluster
                            new_mc1 = temp$new_mc
                            ave_CP_size_mc1[it] = ave_CP_size_mc1[it] + temp$ave_CP_size
                            temp = flipClusterSpins_V2(cluster_mc2, mc2, labelLabel_mc2, rr, cc)
                            cluster_mc2 = temp$cluster
                            new_mc2 = temp$new_mc
                            ave_CP_size_mc2[it] = ave_CP_size_mc2[it] + temp$ave_CP_size
                            
                            mc1 = new_mc1
                            mc2 = new_mc2
                            
                        }
                    }
                    
                } # end of else beta < 1.0
                
            } # else of sw==1
            
            mc1 = new_mc1
            mc2 = new_mc2
            
            
            #             # debug
            #             # show image at the current sweep
            #             par(mfrow=c(2,2))
            #             image(mc1*256)
            #             image(mc2*256)
            
            h1 = 0
            h2 = 0
            for (rr in 1:nrow) {
                for (cc in 1:ncol) {
                    
                    neighbors = findNeighbors(cc, rr, nrow, ncol)
                    
                    # sufficient statistics
                    h1 = h1 + sum(mc1[neighbors] != mc1[rr,cc])
                    h2 = h2 + sum(mc2[neighbors] != mc2[rr,cc])
                }
            }
            # divide by 2 to account for double counting
            # since we are visiting each site twice (via the neighbors)
            h_mc1[it] = 1/(2*n^2) * h1 * 0.5
            h_mc2[it] = 1/(2*n^2) * h2 * 0.5
            
            # check for convergence
            if (abs(h_mc1[it] - h_mc2[it]) < diff_eps || it == maxIter) {
                taus[bb] = it
                
                # 1. plot the sufficient statistics
                df1 = data.frame(iter = rep(1:it, 2), h=c(h_mc1[1:it],h_mc2[1:it]), 
                                 mc=as.factor(rep(1:2,each=it)))
                cat('\n')
                print(ggplot(df1, aes(iter,h, colour=mc)) + geom_line() + xlab("iterations") 
                      + ylab("h(x)") + ggtitle(paste0("Problem 1: ", algo, ": beta = ", beta)))
                
                # 3. get the average sizes of the CPs
                # (number of spins/pixels that are flipped together at each sweep) 
                if (sw == 2) {
                    # the content currently in ave_CP_size_mcX is the 
                    # total # sites flipped per sweep
                    ave_CP_size_mc1[1:it] = ave_CP_size_mc1[1:it]/n^2
                    ave_CP_size_mc2[1:it] = ave_CP_size_mc2[1:it]/n^2
                }
                df4 = data.frame(iter = rep(1:it, 2), ave_CP_size=c(ave_CP_size_mc1[1:it],
                                                                    ave_CP_size_mc2[1:it]), 
                                 mc=as.factor(rep(1:2,each=it)))
                print(ggplot(df4, aes(iter,ave_CP_size, colour=mc)) + geom_line() + 
                          xlab("iterations") + ylab("Average CP size") + 
                          ggtitle(paste0("Problem 3: ", algo, ": beta = ", beta)))
                cat('\n')
                
                # 3. get the average sizes of the CPs for this beta
                beta_CP_size_mc1[bb] = mean(ave_CP_size_mc1[1:it])
                beta_CP_size_mc2[bb] = mean(ave_CP_size_mc2[1:it])
                
                # save the sufficient statistics at convergence
                h_at_convergence[bb] = h_mc1[it]
                
                bb = bb + 1
                
                break
            }
        }
    }
    
    # 2. plot the convergence time in sweep
    gibbs = data.frame(beta=betas[1:3], tau=c(70, 202, 581), 
                       method=as.factor(rep("Gibbs Sampler", 3)))
    df2 = data.frame(beta=betas[1:3], tau=taus[1:3], 
                     method=as.factor(rep(algo, 3)))
    df3 = rbind(df2, gibbs)
    cat('\n')
    print(ggplot(df3, aes(beta, tau, colour=method)) + geom_line() + 
              ggtitle(paste0("Problem 2: ", algo, ": Convergence time in sweeps")))
    cat('\n')
    
    df2 = data.frame(beta=betas, tau=taus, method=as.factor(rep(algo, length(betas))))
    df3 = rbind(df2, gibbs)
    print(kable(df3, caption=paste0("Problem 4: ", algo, ": Convergence time in sweeps")))
    cat('\n')
    
    # 3. plot the average sizes of the CPs for each beta
    df6 = data.frame(beta=rep(betas,2), beta_cp_size=c(beta_CP_size_mc1, beta_CP_size_mc2), 
                     mc=as.factor(rep(1:2,each=length(betas))))
    print(ggplot(df6, aes(beta, beta_cp_size, colour=mc)) + geom_line() + 
              ggtitle(paste0(algo, ": Average size of the connected 
                             components at each sweep")))
    
    # output the sufficient statistics at convergence
    df5 = data.frame(beta=betas, h_star=h_at_convergence)
    print(kable(df5, caption=paste0(algo, ": h(x) at convergence")))
    cat('\n')
    
}
