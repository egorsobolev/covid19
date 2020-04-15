import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

url = "https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/" + \
    "csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_{}_global.csv"

deaths = pd.read_csv(url.format('deaths'), index_col=1)
cases = pd.read_csv(url.format('confirmed'), index_col=1)

def get_country_data(country):
    c = cases.loc[country]
    if c.ndim > 1:
        c = c.sum()
    c = c.iloc[3:]
    c.index=pd.to_datetime(c.index, errors="coerce", format="%m/%d/%y")
    d = deaths.loc[country]
    if d.ndim > 1:
        d = d.sum()
    d = d.iloc[3:]
    d.index=pd.to_datetime(d.index, errors="coerce", format="%m/%d/%y")
    return c, d

def cont(y):
    dy = np.diff(y)
    i0 = dy.size-np.argmin(dy[::-1])
    return i0

def fit(y, i0=0, iN=None):

    dy = np.diff(y)
    if iN is None:
        iN = y.size

    y0 = y[i0]

    dy = dy[i0:iN-1]
    t = np.arange(dy.size)+0.5+i0
    ii = dy!=0
    t = t[ii]
    dy = dy[ii] / np.concatenate([[1], np.diff(t)])

    
    y1 = 0.5*(y[i0+1:iN]+y[i0:iN-1])


    z = np.log(dy/y1[ii])
    #t = np.arange(z.size)+0.5+i0


    Mt = t.mean()
    Mz = z.mean()
    Mtt = (t*t).mean()
    Mzt = (z*t).mean()

    '''
    (pt + q - z) * t = 0
    (pt + q - z) = 0

    p tt + qt = zt
    p t  + q  = z

    D  = <t*t>*<1> - <t>*<t>
    Dp = <z*t>*<1> - <t>*<z>
    Dq = <t*t>*<z> - <zt>*<t>

    exp(p*t)*exp(q) = b*c*exp(-ct)
    c = -p
    b = -exp(q)/p

    '''

    D = Mtt - Mt*Mt
    p = (Mzt - Mt*Mz) / D
    q = (Mtt*Mz - Mzt*Mt) / D

    c = -p
    b = -np.exp(q)/p

    y2 = y[i0:iN]
    t2 = np.arange(y2.size)+i0

    a = np.mean(y2*np.exp(b*np.exp(-c*t2)))
    
    return a,b,c

def fit_all_inv(y, inv):
    prm = []
    for i0,iN in inv:
        a,b,c=fit(y,i0,iN)
        t_max = np.log(b) / c
        r_max = a*np.exp(-1)*c
        f_max = i0<=t_max and t_max <=iN
        print(a,b,c,i0,iN)
        prm.append([a,b,c,i0,iN,t_max,r_max,f_max])
    return prm

def show_inv(y, inv):
    t, r = comp_rate(y)

    fig = plt.figure(figsize=[15,5], tight_layout=True)
    ax = fig.gca()
    ax.semilogy(t, r, '.')

    rmn,rmx = r.min(), r.max()

    for i0, iN in inv:
        plt.plot([i0, i0], [rmn, rmx])

    plt.show()

def comp_rate(y):
    dy = np.diff(y)
    t = np.arange(dy.size) + 0.5
    ii, = np.where(dy != 0)
    r = 2*dy[ii]/(y[ii+1]+y[ii])
    t = t[ii]
    return t, r


def show_model(y, prm, xmax=120):
    fig = plt.figure(figsize=[15,12], tight_layout=True)
    ax = fig.add_subplot(3,1,1)
    ax.set_title('Cumulative cases')
    ax.semilogy(y, '.')
    for cc, [a,b,c,i0,iN,t_max,r_max,f_max] in enumerate(prm):
        t = np.arange(i0,iN+1)
        f = a*np.exp(-b*np.exp(-c*t))
        ax.semilogy(t,f,c=f'C{cc}')
        if f_max:
            ax.semilogy(t_max, a*np.exp(-1),'o',c=f'C{cc}')
            
            
    ax.set_xlim(0,xmax)


    ax = fig.add_subplot(3,1,2)
    ax.set_title('Daily cases')
    dy = np.diff(y)
    t = np.arange(dy.size) + 0.5
    ax.semilogy(t,dy, '.')
    for cc, [a,b,c,i0,iN,t_max,r_max,f_max] in enumerate(prm):
        t = np.arange(i0,iN+1)
        f = a*np.exp(-b*np.exp(-c*t))
        r = b*c*np.exp(-c*t)
        ax.semilogy(t,f*r,c=f'C{cc}')
        if f_max:
            ax.semilogy(t_max, r_max,'o',c=f'C{cc}')

    cc += 1
    t = np.arange(iN+1, xmax)
    f = a*np.exp(-b*np.exp(-c*t))
    r = b*c*np.exp(-c*t)
    ax.semilogy(t, f*r, c=f'C{cc}')
    t_max = np.log(b) / c
    r_max = a*np.exp(-1)*c
    if iN+1<=t_max and t_max <=xmax:
        ax.semilogy(t_max, r_max, 'o', c=f'C{cc}')
    ax.set_xlim(0,xmax)


    ax = fig.add_subplot(3,1,3)
    ax.set_title('Rate')
    t, r = comp_rate(y)
    ax.semilogy(t,r, '.')

    for cc, [a,b,c,i0,iN,t_max,r_max,f_max] in enumerate(prm):
        t = np.arange(i0,iN+1)
        r = b*c*np.exp(-c*t)
        plt.semilogy(t,r, c=f'C{cc}')
    ax.set_xlim(0,xmax)



    plt.show()