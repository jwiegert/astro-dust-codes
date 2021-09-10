# This function just sets sstandard settings for my plots.
#
def settickslabels(xlabel,ylabel,xscale,yscale):
    #
    # Import required packages
    #
    import matplotlib.pyplot as plt
    from matplotlib import rc
    rc('font',**{'family':'serif','serif':['serif']})
    rc('text', usetex=True)
    #
    # Set size of ticks and tick fonts, labels and scales.
    #
    rc('xtick.major',size=8)
    rc('xtick.minor',size=4)
    rc('ytick.major',size=8)
    rc('ytick.minor',size=4)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    plt.xlabel(xlabel,fontsize=18)
    plt.ylabel(ylabel,fontsize=18)
    if xscale == "log":
        plt.scale('log')
    if yscale == "log":
        plt.yscale('log')
