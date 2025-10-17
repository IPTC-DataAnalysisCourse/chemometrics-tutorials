import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import plotly.graph_objs as go
import matplotlib
from matplotlib.colors import Normalize
import matplotlib.cm as cm
from matplotlib.collections import LineCollection
from mpl_toolkits.axes_grid1 import make_axes_locatable

__authors__ = ['gscorreia89', 'flsoares', 'kopeckylukas', 'Hummashazi']
__date__ = "2025/10/17"

def manhattan_plot(pvalues, beta, sig=0.05, instrument='nmr', xaxis=None, yaxis=None):
    """

    :param np.ndarray pvalues: Numpy array with the p-values. These can also be adjusted, or fdr corrected/estimates
    :param np.ndarray beta: Numpy array with the regression coefficients or other effect size estimate.
    :param float sig: Significance value to plot
    :param np.ndarray xvalues: Variable names
    :return: Matplotlib figure with the Manhattan plot
    """
    
    logged_p = -np.log10(pvalues)
    fig, ax = plt.subplots()
    ax.set_title("Manhattan plot")
    
    if instrument == 'nmr':

        ax.set_ylabel(r"Sign($\beta$) $\times$ - $log_{10}$p-value")
        ax.set_xlabel(r"$\delta$ppm")
        if xaxis is None:
            xaxis = np.arange(pvalues.size)
        scatter_plot = ax.scatter(xaxis, np.sign(beta) *logged_p, s=10, c=beta)
        ax.axhline(-np.log10(sig), linestyle='--')
        ax.axhline(- 1*-np.log10(sig), linestyle='--')
    
        # plt.plot(np.mean(X, axis=0).T)
        fig.colorbar(scatter_plot)
        ax.invert_xaxis()
        plt.show()
        
    elif instrument == 'lcms':
        
        ax.set_ylabel('Mass to charge ratio (m/z)')
        ax.set_xlabel('Retention Time')
        
        colormap=plt.cm.coolwarm
        
        # Define colorbar and colormap limits
        maxval = np.max([np.abs(np.max(beta)), np.abs(np.min(beta))])
        maxcol = maxval
        mincol = -maxval
        new_cmap = _shiftedColorMap(colormap, start=0, midpoint=1 - maxcol/(maxcol + np.abs(mincol)), stop=1, name='new')   
        
        # Group significant values
        intensity = np.sign(beta) *logged_p  # The vector for grouping in colormap
        group = np.zeros(np.shape(pvalues))
        group[intensity>-np.log10(sig)] = 1
        group[intensity<- 1*-np.log10(sig)] = 1
        
        # Create a scatter plot with a colormap based on the third vector
        ax = plt.scatter(x=xaxis[group == 0], y=yaxis[group == 0], color='gray', s =10)
        # Show only the selected group of points with a different marker style or color
        ax = sns.scatterplot(x=xaxis[group == 1], y=yaxis[group == 1], hue=beta[group == 1], palette=new_cmap)
        # Customize the color bar
        norm = Normalize(vmin=mincol, vmax=maxcol)
        sm = plt.cm.ScalarMappable(cmap=new_cmap, norm=norm)
        sm.set_array([])
        # Remove the legend and add a colorbar
        ax.get_legend().remove()
        cbar = ax.figure.colorbar(sm)
            
        cbar.set_label(r"Sign($\beta$) $\times$ - $log_{10}$p-value")


def interactive_manhattan(pvalues, beta, sig=0.05, instrument='nmr', xaxis=None, yaxis=None):
    """

    :param np.ndarray pvalues: Numpy array with the p-values. These can also be adjusted, or fdr corrected/estimates
    :param np.ndarray beta: Numpy array with the regression coefficients or other effect size estimate.
    :param float sig: Significance value to plot
    :param np.ndarray xvalues: Variable names
    :return: Data structure ready for interactive plotting with plotly
    """

    data = []
    logged_p = -np.log10(pvalues)

    if instrument == 'nmr':
        if xaxis is None:
            xaxis = np.arange(pvalues.size)
    
        yvals = np.sign(beta) * logged_p
    
        W_str = ["%.4f" % i for i in beta]  # Format text for tooltips
        maxcol = np.max(abs(beta))
    
        Xvals = xaxis
        hovertext = ["ppm: %.4f; W: %s" % i for i in zip(Xvals, W_str)]  # Text for tooltips
    
        point_text = ["p-value: " + pval for pval in pvalues.astype(str)]
    
        manhattan_scatter = go.Scattergl(
            x=Xvals,
            y=yvals,
            mode='markers',
            marker=dict(color=beta, size=5, colorscale='RdBu',
                        cmin=-maxcol, cmax=maxcol, showscale=True),
            text=point_text)
    
        data.append(manhattan_scatter)
    
        xReverse = 'reversed'
        Xlabel = chr(948) + 'ppm 1H'
        Ylabel = r"Sign($\beta$) $\times$ - $log_{10}$p-value"
    
        # Add annotation
        layout = {
            'xaxis': dict(
                title=Xlabel,
                autorange=xReverse),
            'yaxis': dict(title=Ylabel),
            'title': 'Manhattan plot',
            'hovermode': 'closest',
            'bargap': 0,
            'barmode': 'stack',
            'shapes': [{
                'type': 'line',
                'x0': min(Xvals),
                'y0': np.log10(sig),
                'x1': max(Xvals),
                'y1': np.log10(sig),
                'line': {
                    'color': 'rgb(50, 171, 96)',
                    'width': 4,
                    'dash': 'dashdot'}},
                {
                    'type': 'line',
                    'x0': min(Xvals),
                    'y0': -np.log10(sig),
                    'x1': max(Xvals),
                    'y1': -np.log10(sig),
                    'line': {
                        'color': 'rgb(50, 171, 96)',
                        'width': 4,
                        'dash': 'dashdot'}}]}
        fig = {
            'data': data,
            'layout': layout,
        }
    elif instrument == 'lcms':
        
        # Define colorbar and colormap limits
        maxval = np.max([np.abs(np.max(beta)), np.abs(np.min(beta))])
        maxcol = maxval
        mincol = -maxval
        
        # Group significant values
        intensity = np.sign(beta) *logged_p  # The vector for grouping in colormap
        group = np.zeros(np.shape(pvalues))
        group[intensity>-np.log10(sig)] = 1
        group[intensity<- 1*-np.log10(sig)] = 1
    
        # First plot - unselected Features
        manhattan_scatter = go.Scattergl(name="non-Significant",
            x=xaxis[group == 0],
            y=yaxis[group == 0],
            mode='markers',
            opacity=0.5,
            marker=dict(color='gray', size=5))
        
        # Second plot - selected Features
        manhattan_scatter_sel = go.Scattergl(name="Significant",
            x=xaxis[group == 1],
            y=yaxis[group == 1],
            mode='markers',
            opacity=1,
#             marker=dict(color=beta[group == 1], size=10, colorscale=[[0, 'darkblue'], [0.5, 'cornsilk'], [1, 'darkred']],
            marker=dict(color=beta[group == 1], size=10, colorscale='RdBu_r',
                        cmin=-maxcol, cmax=maxcol, showscale=True))

        # Append data
        data = [manhattan_scatter,manhattan_scatter_sel]

        # Create labels
        Xlabel = 'Mass to charge ratio (m/z)'
        Ylabel = r"Sign($\beta$) $\times$ - $log_{10}$p-value"
    
        # Add annotation
        layout = {
            'xaxis': dict(title=Xlabel),
            'yaxis': dict(title=Ylabel),
            'title': 'Manhattan plot',
            'hovermode': 'closest',
            'bargap': 0,
            'barmode': 'stack',
            'legend': dict(
                yanchor="top",
                y=0.99,
                xanchor="right",
                x=0.99)}
        fig = {
            'data': data,
            'layout': layout,
        }

    return fig

def _lineplots(mean, error=None, xaxis=None,color=None,linestyle=None):

    fig, ax = plt.subplots()
    if xaxis is None:
        ax.plot(mean,color=color, linestyle=linestyle)
        xaxis = range(mean.size)
    else:
        ax.plot(xaxis, mean, color=color, linestyle=linestyle)
    if error is not None:
        ax.fill_between(xaxis, mean - error, mean + error, alpha=0.2, color='red')
    return fig, ax


def _barplots(mean, error=None, xaxis=None):

    fig, ax = plt.subplots()
    if xaxis is None:
        xaxis = range(mean.size)

    if error is None:
        ax.bar(xaxis, height=mean)
    else:
        ax.bar(xaxis, height=mean, yerr=error)
    return fig, ax


def _scatterplots(mean, xaxis=None, yaxis=None, colormap=plt.cm.RdYlBu_r, xlabel='Retention Time',
                 ylabel='Mass to charge ratio (m/z)', cbarlabel='Magnitude', marker_size=None, alpha=None):
    """

    """

    colormap = colormap
    maxval = np.max([np.abs(np.max(mean)), np.abs(np.min(mean))])
    maxcol = maxval
    mincol = -maxval
    new_cmap = _shiftedColorMap(colormap, start=0, midpoint=1 - maxcol/(maxcol + np.abs(mincol)), stop=1, name='new')

    fig, ax = plt.subplots()
    # To set the alpha of each point to be associated with the weight of the loading, generate an array where each row corresponds to a feature, the
    # first three columns to the colour of the point, and the last column to the alpha value
    # Return the colours for each feature
    norm = Normalize(vmin=mincol, vmax=maxcol)
    cb = cm.ScalarMappable(norm=norm, cmap=new_cmap)
    cVectAlphas = np.zeros((mean.shape[0], 4))
    cIX = 0
    for c in mean:
        cVectAlphas[cIX, :] = cb.to_rgba(mean[cIX])
        cIX = cIX + 1

    # Set the alpha (min 0.2, max 1)
    cVectAlphas[:, 3] = (((abs(mean) - np.min(abs(mean))) * (1 - 0.2)) / (np.max(abs(mean)) - np.min(abs(mean)))) + 0.2
    if any(cVectAlphas[:, 3] > 1):
        cVectAlphas[cVectAlphas[:, 3] > 1, 3] = 1

    # Plot
    if xaxis is None and yaxis is None:
        raise TypeError("Please, inform xaxis and yaxis")
    elif xaxis is not None and yaxis is None:
        raise TypeError("Please, inform yaxis")
    elif xaxis is None and yaxis is not None:
        raise TypeError("Please, inform xaxis")
    else:   
        ax.scatter(xaxis, yaxis, color=cVectAlphas, s=marker_size, alpha=alpha)
        cb.set_array(mean)
        ax.set_xlim([min(xaxis)-1, max(xaxis)+1])
    
        divider = make_axes_locatable(ax)
        cax = divider.append_axes('right', size='5%', pad=0.05)
        cbar = plt.colorbar(cb, cax=cax)
        cbar.set_label(cbarlabel)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
    return None


def plotLoadings(pcaLoadings, ppm, spectra, title='', figures=None, savePath=None, figureFormat='png', dpi=72, figureSize=(11, 7), cbarlabels ='Loadings'):
    """
    Plot PCA loadings for each component in PCAmodel. For NMR data plots the median spectrum coloured by the loading. For MS data plots an ion map (rt vs. mz) coloured by the loading.
    :param pcaLoadings: Loading vector
    :param ppm: ppm vector to use in x axis
    :param str title: Title for the plot
    :param dict figures: If not ``None``, saves location of each figure for output in html report (see multivariateReport.py)
    """

    Xvals = ppm
    Xlabel = chr(948)+ '1H'

    Yvals = np.median(spectra, axis=0)
    Ylabel = 'Median Intensity'

    cVect = pcaLoadings
    orig_cmap = plt.cm.RdYlBu_r # Red for high, Blue for negative, and we will have a very neutral yellow for 0
    maxval = np.max([np.abs(np.max(cVect)), np.abs(np.min(cVect))])
    maxcol = maxval#npy.max(cVect) # grab the maximum
    mincol = -maxval#npy.min(cVect) # Grab the minimum
    new_cmap = _shiftedColorMap(orig_cmap, start=0, midpoint=1 - maxcol/(maxcol + np.abs(mincol)), stop=1, name='new')

    fig, ax = plt.subplots(figsize=figureSize, dpi=dpi)

    ax.set_rasterized(True)

    lvector = cVect
    points = np.array([Xvals, Yvals]).transpose().reshape(-1,1,2)
    segs = np.concatenate([points[:-1],points[1:]],axis=1)

    cb = LineCollection(segs, cmap=new_cmap)
    cb.set_array(lvector)
    plt.gca().add_collection(cb) # add the collection to the plot
    plt.xlim(Xvals.min()-0.4, Xvals.max() + 0.4) # line collections don't auto-scale the plot
    plt.ylim(Yvals.min()*1.2, Yvals.max()*1.2)
    plt.gca().invert_xaxis()

    cbar = plt.colorbar(cb)
    cbar.set_label(cbarlabels)
    ax.set_xlabel(Xlabel)
    ax.set_ylabel(Ylabel)

    if savePath:
        if figures is not None:
            saveTemp = title + 'PCAloadings'
            figures[saveTemp] = savePath + saveTemp + '.' + figureFormat
        else:
            saveTemp = ''
            plt.savefig(savePath + saveTemp + '.' + figureFormat, bbox_inches='tight', format=figureFormat, dpi=dpi)
            plt.close()

    else:
        plt.show()

    if figures is not None:
        return figures


def _shiftedColorMap(cmap, start=0, midpoint=0.5, stop=1.0, name='shiftedcmap'):
    '''
    From Paul H at Stack Overflow
    http://stackoverflow.com/questions/7404116/defining-the-midpoint-of-a-colormap-in-matplotlib
    Function to offset the "center" of a colormap. Useful for
    data with a negative min and positive max and you want the
    middle of the colormap's dynamic range to be at zero

    Input
    -----
      cmap : The matplotlib colormap to be altered
      start : Offset from lowest point in the colormap's range.
          Defaults to 0.0 (no lower ofset). Should be between
          0.0 and `midpoint`.
      midpoint : The new center of the colormap. Defaults to
          0.5 (no shift). Should be between 0.0 and 1.0. In
          general, this should be  1 - vmax/(vmax + abs(vmin))
          For example if your data range from -15.0 to +5.0 and
          you want the center of the colormap at 0.0, `midpoint`
          should be set to  1 - 5/(5 + 15)) or 0.75
      stop : Offset from highets point in the colormap's range.
          Defaults to 1.0 (no upper ofset). Should be between
          `midpoint` and 1.0.
    '''
    cdict = {
        'red': [],
        'green': [],
        'blue': [],
        'alpha': []
    }

    # regular index to compute the colors
    reg_index = np.linspace(start, stop, 257)

    # shifted index to match the data
    shift_index = np.hstack([
        np.linspace(0.0, midpoint, 128, endpoint=False),
        np.linspace(midpoint, 1.0, 129, endpoint=True)
    ])

    for ri, si in zip(reg_index, shift_index):
        r, g, b, a = cmap(ri)
        cdict['red'].append((si, r, r))
        cdict['green'].append((si, g, g))
        cdict['blue'].append((si, b, b))
        cdict['alpha'].append((si, a, a))

    newcmap = matplotlib.colors.LinearSegmentedColormap(name, cdict)
    # plt.register_cmap(cmap=newcmap)

    return newcmap

