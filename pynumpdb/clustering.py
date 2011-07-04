from sys import stderr
import numpy as np
from scipy.cluster.vq import vq,kmeans,whiten
from scikits.learn import gmm
from _pca import pca_train

def version():
  print """
  pdb_gmm ver. 0.1
  coded by Shiqiao Du in 28,Jun,2010
  kmeans and GMM are supported
  """
def normalize(data,mode="pca",n=10):
  """ normalize and reduce data by PCA"""
  
  if mode == "whiten":
    res = whiten(data)
  elif mode == "pca":
    v,P,res = pca_train(data,n,0,1) 
    print v
    print "eigen ratio is ",v[n-1] / v[0] 
  elif mode == "pca_whiten": 
    v,P,proj = pca_train(data,n,0,1) 
    res = whiten(proj)
  else:
    res = np.array(data)    
 
  return res

# num of variables of gmm
def complexity(data,nc):
  n,ndim = data.shape
  k = nc # num of clust
  k += nc * ndim # num of mean
  k+= nc * (ndim  * (ndim -1) * 0.5 + ndim) # num of cov
  comp = k 
  return comp

# Akaike information criterion
def AIC(log_like,data,nc):
  return (-log_like + complexity(data,nc)) * 2.0

# Bayesian information criterion
def BIC(log_like,data,nc):
  return -log_like * 2.0 + complexity(data,nc) * np.log(len(data))

# Exception for specifying wrong clustering mode
class ClusteringModeError(Exception):
  def __init__(self,mode):
    self.mode = mode
  
  def __str__(self):
    m = "undefined mode '$s'" % self.mode  
    return m

# Main class
class ClusteringMachine:
  
  def __init__(self,nclust):
    self.nclust = nclust # num of clusters
    self.sort_order = False # how to sort members of each cluster 
    self.gm = None # gmm instance

  def clustering(self,data,mode="kmeans",itermax=100,threshold=1.0e-2):
    """
    clustering data
    
    input 
      data : 2d-array, data to be clusterd
      mode : strings, clustering method, "kmeans" or "gmm"
      itermax : int, max cycle of iteration in gmm train
      threshold : float, gmm iteration is converged if changes are below it
      
    return
      codes : array, codes[i] is the corresponding cluster of data[i]
      score : array, in kmeans, score is the distance, in gmm probability
      iter_log : array, log of iteration (only gmm)
      aic : float, AIC of the model (only gmm)
      bic : float, BIC of the model (only gmm)
    """
    
    print "running %s mode" % mode
    
    if mode == "kmeans": # run kmeans mode
      self.sort_order = False # ascending order
      
      # get centroids
      centroids, var = kmeans(data,self.nclust)
      
      # assign data to the nearest cluster
      codes, distances = vq(data,centroids)
      
      return codes,distances
    
    elif mode == "gmm": # run gmm mode
      
      self.sort_order = True # descending order
      
      # construct a GMM instance
      self.gm = gmm.GMM(self.nclust,cvtype="full")
      
      # train gmm 
      iter_log = self.gm.fit(data,itermax,thresh=threshold)
      if len(iter_log) == itermax :
        stderr.write("warning!! EM step not converged\n")
      
      # assign data to the nearest cluster
      logprobs, codes = self.gm.decode(data)
      
      # calc probability that each datum belongs to corresponding cluster 
      lpr = gmm.lmvnpdf(data,self.gm.means, \
                        self.gm._covars,self.gm._cvtype)
      probs = np.array([lpr[i,codes[i]] for i in range(len(lpr))])
      
      # calc AIC and BIC for model evaluation
      aic = AIC(iter_log[-1],data,self.nclust)
      bic = BIC(iter_log[-1],data,self.nclust)
      
      return codes,probs,iter_log,aic,bic
    
    else:
      raise ClusteringModeError(mode)

  def rankCodes(self,codes):
    """
    rank codes according to its population
    
    input 
      codes : array

    return
      new_codes : array, renumbered codes
      pop : array, population of sorted codes
    """
    
    # sort clusters by their size
    icount = zip(range(self.nclust),np.bincount(codes))
    icount.sort(reverse=True,key = lambda x:int(x[1]))

    # calc population of each cluster
    totn = 1.0/sum(ic[1] for ic in icount)
    pop = [ic[1] * totn for ic in icount]

    # reassign each datum to the sorted clusters
    rank = {}
    for i,r in enumerate(icount):
      rank[r[0]] = i
    new_codes = [rank[codes[i]] for i in range(len(codes))]

    return np.array(new_codes), np.array(pop)

  def sortedCluster(self,codes,scores,order=0):
    """ 
    assing indecies of data to cluster

    input
      codes : array
      scores : array, sort members of each cluster 
              according to this value
      order : int (optional), -1: increscentm, 1 decrescent order, 0 auto
      
    return
      clusters : array of cluster, members are sorted by scores
    """

    # set sort order if needed
    if order > 0:
      self.sort_order = True # larger is better
    elif order < 0:
      self.sort_order = False # smaller is better
    else :
      pass

    # make sets of i and scores[i]
    clusters = [[] for i in range(self.nclust)]
    for i in range(len(codes)):
      clusters[codes[i]].append((i,scores[i]))
    
    # sort by score
    for c in clusters:
      c.sort(key=lambda x:float(x[1]),reverse=self.sort_order)
      
    return clusters

