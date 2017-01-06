import warnings
import numpy as np
import pandas as pd
import scipy.stats as st
import statsmodels as sm
import matplotlib
import matplotlib.pyplot as plt
import pickle
import time

matplotlib.rcParams['figure.figsize'] = (16.0, 12.0)
matplotlib.style.use('classic')

# Create models from data
def best_fit_distribution(data, bins=200, ax=None):
    """Model data by finding best fit distribution to data"""
    # Get histogram of original data
    count, xbin_edges = np.histogram(data, bins=bins, normed=True)

    # Distributions to check
    DISTRIBUTIONS = [        
        st.alpha,st.anglit,st.arcsine,st.beta,st.betaprime,st.bradford,st.burr,st.cauchy,st.chi,st.chi2,st.cosine,
        st.dgamma,st.dweibull,st.erlang,st.expon,st.exponnorm,st.exponweib,st.exponpow,st.f,st.fatiguelife,st.fisk,
        st.foldcauchy,st.foldnorm,st.frechet_r,st.frechet_l,st.genlogistic,st.genpareto,st.gennorm,st.genexpon,
        st.genextreme,st.gausshyper,st.gamma,st.gengamma,st.genhalflogistic,st.gilbrat,st.gompertz,st.gumbel_r,
        st.gumbel_l,st.halfcauchy,st.halflogistic,st.halfnorm,st.halfgennorm,st.hypsecant,st.invgamma,st.invgauss,
        st.invweibull,st.johnsonsb,st.johnsonsu,st.ksone,st.kstwobign,st.laplace,st.levy,st.levy_l,st.levy_stable,
        st.logistic,st.loggamma,st.loglaplace,st.lognorm,st.lomax,st.maxwell,st.mielke,st.nakagami,st.ncx2,st.ncf,
        st.nct,st.norm,st.pareto,st.pearson3,st.powerlaw,st.powerlognorm,st.powernorm,st.rdist,st.reciprocal,
        st.rayleigh,st.rice,st.recipinvgauss,st.semicircular,st.t,st.triang,st.truncexpon,st.truncnorm,st.tukeylambda,
        st.uniform,st.vonmises,st.vonmises_line,st.wald,st.weibull_min,st.weibull_max,st.wrapcauchy
    ]

    # Best holders
    best_distribution = st.norm
    best_params = (0.0, 1.0)
    best_sse = np.inf
    StartTime = time.time()
    # Estimate distribution parameters from data
    for i in range(len(DISTRIBUTIONS)):

        # Try to fit the distribution
        try:
            # Ignore warnings from data that can't be fit
            with warnings.catch_warnings():
                warnings.filterwarnings('ignore')

                # fit dist to data
                params = distribution[i].fit(data, floc = 0)

                # Separate parts of parameters
                arg = params[:-2]
                loc = params[-2]
                scale = params[-1]

                # Calculate fitted PDF and error with fit in distribution
                pdf = distribution[i].pdf(xbin_edges, *params)
                sse = np.sum(np.power(count - pdf, 2.0))

                # identify if this distribution is better
                if best_sse > sse > 0:
                    best_distribution = distribution[i]
                    best_params = params
                    best_sse = sse
        except Exception:
            pass
        statusbar(i,len(DISTRIBUTIONS),StartTime=StartTime,Title='best_fit_distribution')

    return (best_distribution.name, best_params)

def make_pdf(dist, params, size=10000):
    """Generate distributions's Propbability Distribution Function """

    # Separate parts of parameters
    arg = params[:-2]
    loc = params[-2]
    scale = params[-1]

    # Get sane start and end points of distribution
    start = dist.ppf(0.01, *arg, loc=loc, scale=scale) if arg else dist.ppf(0.01, loc=loc, scale=scale)
    end = dist.ppf(0.99, *arg, loc=loc, scale=scale) if arg else dist.ppf(0.99, loc=loc, scale=scale)

    # Build PDF and turn into pandas Series
    x = np.linspace(start, end, size)
    y = dist.pdf(x, loc=loc, scale=scale, *arg)
    pdf = pd.Series(y, x)

    return pdf

def statusbar(i,N,**kwargs):
    """
    i is the current iteration (must be an int) and N is the length of 
    the range (must be an int). i must also be in [0,N). 
    
    ~~~~~~~~~~~~~~
    **kwargs
    ~~~~~~~~~~~~~~
    
    StartTime should equal time.time() and should be defined before your
    loop to ensure that you get an accurate representation of elapsed time.

    Title should be a str that will be displayed before the statusbar. Title
    should be no longer than 25 characters.

    ~~~~~~~~~~~~~~

    NOTE: you should place a print('\n') after the loop to ensure you
    begin printing on the next line.

    """
    import time
    StartTime = kwargs.get("StartTime",False)
    Title = kwargs.get("Title",'')

    assert type(i)==int, "i must be an int"
    assert type(N)==int, "N must be an int"
    assert N>i, "N must be greater than i"
    assert N>0, "N must be a positive integer"
    assert i>=0, "i must not be negative (can be zero)"
    assert type(Title) == str, "Title should be a string"
    assert len(Title) <= 25, "Title should be less than 25 characters"
    if Title != '': Title = ' '*(25-len(Title)) + Title + ': '
    statusbar = Title +'[' + '\u25a0'*int((i+1)/(N/50)) + '\u25a1'*(50-int((i+1)/(N/50))) + '] '
    if StartTime != False:
        print(statusbar + '{0:1.1f}'.format((i+1)/N*100) + '% complete, ' + '{0:1.1f}'.format(time.time() - StartTime) + 'sec        \r', end='')
    else:
        print(statusbar + '{0:1.1f}'.format((i+1)/N*100) + '% complete           \r',end = '')
# Load data from statsmodels datasets
data = pickle.load(open("2804ErrorDistribution.pkl","rb"))
data.remove(0.0)
data = np.array(data)
bins = 50
# # Plot for comparison
# plt.figure(figsize=(12,8))
# plt.hist(data/5501,normed = 1, bins = bins)
# ax = plt.gca()
# # Save plot limits
# dataYLim = ax.get_ylim()

# plt.figure(figsize=(12,8))
# plt.hist(data/5501,normed = 1, bins = bins,histtype='step',cumulative=True,color = 'k',lw = 1)
# ax = plt.gca()
# # Save plot limits
# dataYLim = ax.get_ylim()

# # Find best fit distribution
# best_fit_name, best_fir_paramms = best_fit_distribution(data, bins, ax)
# best_dist = getattr(st, best_fit_name)

# # Update plots
# ax.set_ylim(dataYLim)
# ax.set_title('Error Distribution For Trial 2804')
# ax.set_xlabel('Sum of Error')
# ax.set_ylabel('Frequency')

# # Make PDF
# pdf = make_pdf(best_dist, best_fir_paramms)

# # Display
# plt.figure(figsize=(12,8))
# ax = pdf.plot(lw=2, label='PDF', legend=True)
# plt.hist(data,normed = True, bins=bins)

# param_names = (best_dist.shapes + ', loc, scale').split(', ') if best_dist.shapes else ['loc', 'scale']
# param_str = ', '.join(['{}={:0.2f}'.format(k,v) for k,v in zip(param_names, best_fir_paramms)])
# dist_str = '{}({})'.format(best_fit_name, param_str)

# ax.set_title(u'Error Distribution For Trial 2804\nwith best fit distribution \n' + dist_str)
# ax.set_xlabel('Sum of Error')
# ax.set_ylabel('Frequency')

plt.figure(figsize=(12,8))
gamma = st.gamma
a,loc,b = gamma.fit(data/5501)
param = (a,loc,b)
x = np.linspace(0,(data/5501).max(),10000)
pdf_fitted = gamma.pdf(x,*param)
plt.plot(x,pdf_fitted,'r',lw=2)
plt.hist(data/5501,normed=True,bins=50,facecolor = '#BFD6D5')
ax1 = plt.gca()
ax2 = ax1.twinx()
cdf_fitted = gamma.cdf(x,param[0],loc=param[1],scale=param[2])
ax2.plot(x,cdf_fitted,'k',lw=3)
plt.show()