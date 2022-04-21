import numpy as np
from datetime import datetime
import time as time
import matplotlib.pyplot as plt
import scipy.interpolate as spi

def imputeMissing(X, method, parameters=None):
    """ Impute missing features of X (NaNs) in one of several simple ways
    X,parameters = imputeMissing(X, method, parameters) 
    Missing values are denoted by NaN
    methods are:
      'constant' : fill with a constant value
      'mean'     : fill all missing values with the mean over that feature
      'median'   : fill "" with the median value
      'gaussian' : fill with the conditional mean assuming a Gaussian model on X (w/ shrinkage to N(0,1))
      'gausscopula' : 'gaussian' but with univariate copulas to attempt to preserve marginal distributions
    parameters : (optional) method-specific information to use in imputation:
      'constant' : the constant value to fill with
      'mean', 'median' : a vector of values (one per feature) to fill with
      'gaussian' : (mean,Covar), the mean and covariance to use for the Gaussian
      'gausscopula' : (mean,Covar,Empiricals), also includes empirical marginal CDFs
    """
    X = X.copy()
    m,n = X.shape
    method = method.lower()
    def nanEval(X, lam):
       e = np.zeros( (X.shape[1],) )
       for i in range(X.shape[1]):
         e[i] = lam(X[ ~np.isnan(X[:,i]),i ])
       return e
    # First, create imputation parameters if not provided:
    if parameters is None:
        if method == 'mean':
            #parameters = np.nanmean(X, axis=0)
            parameters = nanEval(X, lambda X: np.mean(X))
        if method == 'median':
            #fillValue = np.nanmedian(X, axis=0)
            parameters = nanEval(X, lambda X: np.median(X))
        if method == 'gaussian':
            mu = nanEval(X, lambda X: np.mean(X))
            for i in range(n):
                mi = float( np.sum(~np.isnan(X[:,i])) )
                mu[i] *= mi/(mi+n)  # shrink mean toward zero by m counts
            cov = np.zeros((n,n))
            for i in range(n):
                for j in range(i,n):
                    nans = np.isnan(X[:,i]) | np.isnan(X[:,j])
                    mij  = float( np.sum(~nans) )
                    cov[i,j] = np.mean( (X[~nans,i]-mu[i])*(X[~nans,j]-mu[j]) )
                    # cov[i,j] *= mij/(mij+n)         # shrink towards
                    cov[i, j] *= mij / (mij + 5)  # shrink towards
                    # if i==j: cov[i,j] += n/(mij+n)  #  identity matrix
                    if i==j: cov[i,j] += 5/(mij+5)  #  identity matrix
                    cov[j,i] = cov[i,j]
            parameters = mu,cov
    # Now, apply imputation paramters to fill in the missing values
    if method == 'constant':
        X[ np.isnan(X) ] = parameters
    if method == 'mean' or method == 'median':
        for i in range(n):
            X[ np.isnan(X[:,i]), i] = parameters[i]
    if method == 'gaussian':
        mu,Sig = parameters
        for j in range(m):
            nans = np.argwhere(np.isnan(X[j,:])).flatten()
            oks  = np.argwhere(~np.isnan(X[j,:])).flatten()
            X[j,nans] = mu[nans] + Sig[np.ix_(nans,oks)].dot( np.linalg.inv(Sig[np.ix_(oks,oks)]).dot( (X[j,oks]-mu[oks]).T ) ).T
    if method == 'gausscopula':  # TODO: don't ignore input parameters?
        X,PEmp = toGaussCopula(X);
        X,params = imputeMissing(X,'gaussian');
        parameters = params[0],params[1],PEmp
        X = fromGaussCopula(X,PEmp);
    return X,parameters

def toGaussCopula(X):
    '''Convert data X into Gaussian-like marginals via a copula transformation'''
    import scipy.stats
    XC = np.copy(X)
    PEmp = []
    for c in range(X.shape[1]):
        mask = ~np.isnan( X[:,c] ); m = np.sum(mask);
        srt  = np.argsort(X[mask,c]); rank = 0*srt.astype(float);
        rank[srt] = np.linspace(1./(m+1),1.*m/(m+1),m); #(+1.)*1./(np.sum(mask)+1)
        PEmp.append( np.sort(X[mask,c]) )
        XC[mask,c] = scipy.stats.norm.ppf(rank)
    return XC,PEmp

def fromGaussCopula(XC,PEmp):
    """Convert from standard Gaussian marginals to empirical marginals via copula xform"""
    import scipy.stats
    X = np.copy(XC)
    for c in range(XC.shape[1]):
        rank = scipy.stats.norm.cdf(XC[:,c])
        m = len(PEmp[c])
        X[:,c] = np.interp( rank, np.linspace(1./(m+1),1.*m/(m+1),m), PEmp[c] )
    return X

