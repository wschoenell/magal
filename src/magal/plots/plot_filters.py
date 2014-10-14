import numpy as np
import matplotlib.pyplot as plt
from magal.io.readfilterset import FilterSet

def plot_filter_transmissions(filename, filterset, ccd, exclude_filters=(None,), sdss_path=None, ax=None, ylabel=True):
    filter_file = filename

    if not ax:
        plt.figure(1)
        plt.clf()
        ax = plt.axes()

    if sdss_path:
        Filter_sdss = FilterSet(sdss_path)
        Filter_sdss.load('sdss', '1')
        f_sdss = Filter_sdss.filterset
        for i in np.unique(f_sdss['ID_filter']):
            if not i in ['u']:
                ax.fill(f_sdss[f_sdss['ID_filter'] == i]['wl'], f_sdss[f_sdss['ID_filter'] == i]['transm']*100, alpha=.6)

    # Read FilterSet file
    Filter = FilterSet(filename)
    Filter.load(filterset, ccd)
    F = Filter.filterset
    #Define a mean wl vector to the filterset
    print 'Filters: ', np.unique(F['ID_filter'])
    for i in np.unique(F['ID_filter']):
        if not i in exclude_filters:
            ax.plot(F[F['ID_filter'] == i]['wl'], F[F['ID_filter'] == i]['transm']*100) #/F[F['ID_filter'] == i]['transm'].sum())



    ax.set_xlabel('$\lambda \ [\AA]$')
    if ylabel:
        ax.set_ylabel('$R_\lambda\  [\%]$')
    ax.set_xlim(3000, 10000)
    plt.ylim(0, 100)
    ylow = .18
    yupp = plt.ylim()[1]
    ax.set_ylim(ylow, yupp)
    # plt.savefig(thesis_fig_dir+'figure1.eps', format='eps')
