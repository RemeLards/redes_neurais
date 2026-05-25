# https://pandas.pydata.org/pandas-docs/stable/index.html
# from tabulate import tabulate
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
# https://scikit-learn.org/stable/modules/generated/sklearn.preprocessing.LabelBinarizer.html
from sklearn.preprocessing import LabelBinarizer
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import confusion_matrix, classification_report,\
    accuracy_score, f1_score
import numpy as np
import sys
import os
import copy

'''
sys.path.append('/home/thomas/Nextcloud2/book/soft/lib/')
sys.path.append('E:/Nextcloud2/book/soft/lib/')
from LinearMachine import LinearMachine
'''

from sklearn.neighbors import KNeighborsClassifier

def tex_setup(usetex=True):
    if not usetex:
        print('Not setting any Matplotlib parameters for LaTeX ...')
        return
    # import matplotlib.style
    # plt.style.use('classic')
    # plt.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
    # plt.rc('font',**{'family':'serif','serif':['Computer Modern Roman']})
    # for Palatino and other serif fonts use:
    # plt.rc('font',**{'family':'serif','serif':['Palatino']})
    plt.rc('text', usetex=usetex)
    plt.rcParams['text.latex.preamble'] = r'\usepackage{amsmath}'
    plt.rc('font', **{'family': 'serif', 'serif': ['DejaVu Sans']})

    # Set the font size. Either an relative value of
    # 'xx-small', 'x-small', 'small', 'medium', 'large', 'x-large',
    # 'xx-large' or an absolute font size, e.g., 12.
    # https://matplotlib.org/api/font_manager_api.html#matplotlib.font_manager.FontProperties.set_size
    texfigparams = {
        'axes.labelsize': 9,
        # 'axes.linewidth': 0.01,
        'axes.titlesize': 7,
        # 'figure.figsize': (8, 9),
        'font.size': 9,
        'legend.fontsize': 7,
        'legend.loc': 'lower center',
        # 'legend.loc': 'upper left',
        'lines.linewidth': 1,
        # 'savefig.dpi': 600,
        # 'savefig.bbox': 'tight',
        # 'ps.usedistiller': 'ghostscript',
        # 'ps.usedistiller': 'xpdf',
        'text.latex.preamble': r'\usepackage{amsmath}',
        'text.latex.preamble': r'\usepackage{amssymb}',
        'text.latex.preamble': r'\usepackage{bm}',
        # 'text.latex.preamble': r'\usepackage{xcolor}',
        'pgf.preamble': r'\usepackage{amsmath}',
        'pgf.preamble': r'\usepackage{amssymb}',
        'pgf.preamble': r'\usepackage{bm}',
        # 'pgf.preamble': r'\usepackage{xcolor}',
        'xtick.labelsize': 7,
        'ytick.labelsize': 7
    }

    plt.rcParams.update(texfigparams)
    plt.rc('text', usetex=usetex)
    plt.rc('font', family='sans-serif')
    # print('Matplotlib: rcParams=\n', plt.rcParams, file=sys.stdout) ; quit()

# tex_setup()



def read_X(datafn, sep=';'):
    """Read data matrix into a pandas dataframe."""
    df = pd.read_csv(datafn, sep=';')

    # print(df)
    # print('index: ', df.index)
    # print('columns: ', df.columns)
    # print('info: ', df.info)
    # print('shape: ', df.shape); input('...')

    # print(tabulate(df, headers='keys', tablefmt='psql'))
    return df


def dataframe2sklearn(df):
    """Convert a pandas dataframe into sklearn readable format X, y, Y."""
    featname = list(df.columns.values)
    data = df.to_numpy()
    X = data[:, 0:-1]
    labels = data[:, -1]

    # print('data=\n', data, 'shape=', data.shape)
    # print('X=\n', X, 'shape=', X.shape)
    # print('featname=\n', featname)
    # print('labels=\n', labels, 'shape=', labels.shape)
    return X, labels, featname


def ploty_signals(X, labels, featname, mask=None,
                  NORMVAL=np.array([20.0, 0.25, 30.0, 2.00, 2.85, 17.11,
                                    80.0, 0.9, 0.25, 20.0, 56250.0, 25.3, 40.7,
                                    0.9, 0.0, 0.0, 0.0, 0.0])):
    """Plot multichannel signal."""

    def filter_featname(featname, mask):
        if mask is None:
            return featname
        fn = [featname[i] for i in mask]
        return fn

    def filter_vars(X, featname, mask):
        return (copy.copy(X[:, np.array(mask, dtype=int)]),
                filter_featname(featname, mask))

    # plot the multi-channel signal
    n, d = X.shape
    x = np.linspace(0, n-1, n)

    if mask is not None:
        X, featname = filter_vars(X, featname, mask)
        NORMVAL = filter_featname(NORMVAL, mask)
    # print('X.shape=', X.shape, 'featname=', featname)

    n, numsensors = X.shape

    # tex_setup(usetex=usetex)
    fig, ax = plt.subplots(numsensors, 1, sharex=True)

    def cm2inch(cm): return cm/2.54
    widthcm = 24  # Real sizes later in the LaTeX file
    heigthcm = 25
    fig.set_size_inches([cm2inch(widthcm), cm2inch(heigthcm)])
    ax[numsensors-1].set_xlabel('$t [min]$')
    allsignal = X

    # https://matplotlib.org/3.1.1/api/cm_api.html#matplotlib.cm.get_cmap
    colors = plt.colormaps['tab20']

    # Analyze transitions between labels
    condition_change = []
    t = 0
    while t < n-1:
        transition = labels[t] is not labels[t+1]
        if transition:
            # ax.axvline(x=t)
            condition_change.append(t)
            print('t=%5d label transition: %s ===>%s' %
                  (t+1, labels[t], labels[t+1]))
        t += 1
    if condition_change != []:
        print('condition_change=', condition_change)

    for j in range(numsensors):
        # s = numsensors - j - 1
        s = j

        ax[s].plot(x, allsignal[:, s], linewidth=0.5,
                   color=colors(s/numsensors))
        ax[s].axhline(y=NORMVAL[j], linestyle='--', lw=0.5,
                      label=None, color='green')

        # ylabel = '$x_{' + str(s+1) + '} =$' + featname[s]
        # print('featname[%d]=' % s,  featname[s], 'last=', featname[s][-1])
        if featname[s][-1] == '$':
            ylabel = featname[s]
        else:
            ylabel = '$'+featname[s]+'$'
        usetex = True
        ax[s].tick_params(axis='both', which='major', labelsize=7)
        ax[s].tick_params(axis='both', which='minor', labelsize=7)

        ax[s].set_ylabel(ylabel, rotation=0, labelpad=0,
                         verticalalignment='center',
                         horizontalalignment='right',
                         usetex=usetex, fontsize=8)
        for i in range(len(condition_change)):
            xpos = condition_change[i]
            ax[s].axvline(x=xpos, linestyle=':', color='k',
                          linewidth=1.0, label='Transição')
            if s is numsensors-1:
                ylim = ax[s].get_ylim()
                # yoff = ylim[0] - 180.0
                yoff = ylim[0] - 3
                # print('ylim=', ylim)
                transtxt = '$'+str(labels[condition_change[i]+1])+'$'
                # print('s=', s, 's is numsensors-1=', s is numsensors-1,
                #       'transtxt=', transtxt, 'xpos=', xpos, 'yoff=', yoff,
                #       'condition_change[i]=', condition_change[i])
                ax[s].text(xpos, yoff, transtxt,
                           horizontalalignment='center',
                           verticalalignment='top', color='red',
                           rotation=90, fontsize=7)

    # plt.tight_layout()
    # dropfigdir = None
    figname = 'cstr_signals.png'  # .pgf
    print('Saving individual signal plot in ', figname)
    fig.savefig(figname, bbox_inches='tight')
    plt.show()


def pair_plot(df, featname):
    """Plot pairwise 2-D feature space."""
    # Seaborn does not work with np arrays. We need a Dataframe
    df = df.set_axis(featname, axis='columns')

    # print('df.columns=\n', df.columns)

    print('Generating pair plot ...')
    # kde = Kernel Density Estimation, based on Gaussian kernels
    # kernel_wid_1D = 0.5
    # kernel_wid_2D = 1.0
    # https://seaborn.pydata.org/generated/seaborn.pairplot.html
    pp = sns.pairplot(df, hue=featname[-1], diag_kind='kde')  # ,
                      #  diag_kws={'bw_adjust': kernel_wid_1D})
    pp.map_upper(sns.scatterplot)
    pp.map_lower(sns.kdeplot, linewidths=0.5)
    pp.map_diag(sns.kdeplot, linewidth=0.5)
    # pp.map_lower(sns.kdeplot, bw_adjust=kernel_wid_2D, linewidths=0.5)
    # pp.map_diag(sns.kdeplot, bw_adjust=kernel_wid_1D, linewidth=0.5)

    # pp.map_diag(sns.kdeplot, linewidth=0.5)
    figname = 'cstr_pairplot.svg'  # .pgf
    print('Saving pair plot in ', figname)
    plt.savefig(figname, bbox_inches='tight')
    plt.show()
    # pair_plor documentation -> https://seaborn.pydata.org/generated/seaborn.pairplot.html
    # https://seaborn.pydata.org/generated/seaborn.kdeplot.html


def numeric_labels(y):
    """Return numeric labels."""
    # One-Hot if more than two classes, binary for two labels
    classes = np.unique(y)
    # print('numeric_labels> classes=', classes); input('...')
    numclasses = len(classes)
    if numclasses == 2:
        neg_label = -1
    else:
        neg_label = 0
    lb = LabelBinarizer(neg_label=neg_label)
    Y = lb.fit_transform(y)
    # print('Y=\n', Y, 'shape=', Y.shape)
    return Y


def resubstitution(X, y, featname):
    """Classification with resubstitution, Linear Machine classifier."""

    classes = np.unique(y)
    scaler = StandardScaler()
    X = scaler.fit_transform(X)

    '''
    model = LinearMachine()
    model.fit(X, y)
    params = model.get_params(deep=True)
    W, bias, weights = params['W'], params['bias'], params['weights']
    print('W=\n', W)
    print('bias=\n', bias)
    print('weights=\n', weights)
    '''
    model = KNeighborsClassifier(n_neighbors=3)
    model.fit(X, y)
    y_pred = model.predict(X)

    print('\n==========> Resubstitution of training data:\n')
    print('Classification Report for all features and all classes: ')
    print(classification_report(y, y_pred, target_names=classes, digits=3))
    print('Accuracy=', '%.2f %%' % (100*accuracy_score(y, y_pred)))
    print('Confusion Matrix: ')
    print(confusion_matrix(y, y_pred))


def k_fold(X, y, featname, K=10):
    """Classification with k-fold."""

    classes = np.unique(y)
    scaler = StandardScaler()
    skf = StratifiedKFold(n_splits=K, shuffle=False)
    y_pred_overall = []
    y_test_overall = []

    # model = LinearMachine()
    model = KNeighborsClassifier(n_neighbors=3)

    for train_index, test_index in skf.split(X, y):
        # print("TRAIN:", train_index, "TEST:", test_index)
        X_train, X_test = X[train_index], X[test_index]
        y_train, y_test = y[train_index], y[test_index]
        X_train = scaler.fit_transform(X_train)
        X_test = scaler.transform(X_test)

        model.fit(X_train, y_train)
        y_pred = model.predict(X_test)
        y_pred_overall = np.concatenate([y_pred_overall, y_pred])
        y_test_overall = np.concatenate([y_test_overall, y_test])

    accuracy = 100*accuracy_score(y_test_overall, y_pred_overall)
    print('\n==========> K-fold cross validation:\n')
    print('Model ', K, '- Fold Classification Report: ')
    print(classification_report(y_test_overall, y_pred_overall,
                                target_names=classes, digits=3))
    print('Accuracy=', '%.2f %%' % accuracy)
    print('Macro-averaged F1=', '%.3f'
          % (f1_score(y_test_overall, y_pred_overall, average='macro')))
    print('Micro-averaged F1=', '%.3f' %
          (f1_score(y_test_overall, y_pred_overall, average='micro')))
    print('Model Confusion Matrix: ')
    print(confusion_matrix(y_test_overall, y_pred_overall))


def main(datafn=None, mask=None):
    """Execute main program."""

    if datafn is None:
        argc = len(sys.argv)
        # print('Number of arguments:', argc, 'arguments.')
        # print('Argument List:', str(sys.argv))
        if argc > 1:
            datafn = sys.argv[1]
        else:
            print('Usage: %s <datafilename.csv>' % sys.argv[0])
            return

    df = read_X(datafn=datafn)
    X, y, featname = dataframe2sklearn(df)
    ploty_signals(X, y, featname, mask=mask)
    # pair_plot(df, featname)
    resubstitution(X, y, featname)
    # k_fold(X, y, featname)


if __name__ == "__main__":
    datafn = 'X_py.csv' if 'SPYDER_DEV' in os.environ else None
    mask = [0, 1, 9, 10, 13, 14, 17]
    mask = None
    main(datafn=datafn, mask=mask)
