import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
from numpy.linalg import *

def plot1(rho, delta):
    f, axarr = plt.subplots(1, 3)
    axarr[0].set_title('DECISION GRAPH')
    axarr[0].scatter(rho, delta, alpha=0.6, c='red')
    axarr[0].set_xlabel(r'$\rho$')
    axarr[0].set_ylabel(r'$\delta$')
    axarr[1].set_title('DECISION GRAPH 2')
    axarr[1].scatter(np.arange(len(rho)) + 1, -np.sort(-rho * delta), alpha=0.6, c='red')
    axarr[1].set_xlabel('Sorted Sample')
    axarr[1].set_ylabel(r'$\rho*\delta$')
    return (f, axarr)


def plot2(axarr, rho, delta, cmap, cl, icl, XY, NCLUST):
    axarr[0].scatter(rho, delta, alpha=1, c='red')
    axarr[1].scatter(np.arange(len(rho)) + 1, -np.sort(-rho * delta), alpha=1, c='red')
    for i in range(NCLUST):
        axarr[0].scatter(rho[icl[i]], delta[icl[i]], alpha=0.8, c=cmap[i])
        axarr[1].scatter(i + 1, -np.sort(-rho[icl[i]] * delta[icl[i]]), alpha=0.8, c=cmap[i])
    axarr[2].set_title('2D multidimensional scaling')
    axarr[2].scatter(XY[:, 0], XY[:, 1], alpha=0.8, c=cmap[list(cl)])
    axarr[2].set_xlabel('X')
    axarr[2].set_ylabel('Y')


def DCplot(dist, XY, ND, rho, delta, ordrho, dc, nneigh, rhomin, deltamin):
    f, axarr = plot1(rho, delta)

    # global nneigh, ordrho
    #
    # maxd = dist.max()
    # ordrho = (-rho).argsort()
    # delta = np.zeros(ND)
    # nneigh = np.zeros(ND)
    # delta[ordrho[0]] = -1
    # nneigh[ordrho[0]] = 0
    #
    # for ii in range(1, ND):
    #     delta[ordrho[ii]] = maxd
    #     for jj in range(ii):
    #         if dist[ordrho[ii], ordrho[jj]] < delta[ordrho[ii]]:
    #             delta[ordrho[ii]] = dist[ordrho[ii], ordrho[jj]]
    #             nneigh[ordrho[ii]] = ordrho[jj]

    def onclick(event):
        global rhomin, deltamin
        if event.xdata != None and event.ydata != None:
            rhomin = event.xdata
            deltamin = event.ydata
            print('Cutoff: (min_rho, min_delta): (%.2f, %.2f)' % (rhomin, deltamin))
            NCLUST = 0
            cl=[-1]*ND
            # 1000 is the max number of clusters
            icl=[-1]*1000
            for i in range(ND):
                if rho[i] > rhomin and delta[i] > deltamin:
                    cl[i] = NCLUST
                    icl[NCLUST] = i
                    NCLUST = NCLUST + 1

            print('NUMBER OF CLUSTERS: %i' % (NCLUST))
            print('Performing assignation')
            # assignation

            # nneigh1 = nneigh.astype(int)
            # nneigh = nneigh1
            # ordrho1 = ordrho.astype(int)
            # ordrho = ordrho1

            for i in range(ND):
                if cl[ordrho[i]] == -1:
                    cl[ordrho[i]] = cl[nneigh[ordrho[i]]]

            # halo
            # cluster id start from 1, not 0
            ## deep copy, not just reference
            halo = np.zeros(ND)
            halo[:] = cl

            if NCLUST > 1:
                bord_rho = np.zeros(NCLUST)
                for i in range(ND - 1):
                    for j in range((i + 1), ND):
                        if cl[i] != cl[j] and dist[i, j] <= dc:
                            rho_aver = (rho[i] + rho[j]) / 2
                            if rho_aver > bord_rho[cl[i]]:
                                bord_rho[cl[i]] = rho_aver
                            if rho_aver > bord_rho[cl[j]]:
                                bord_rho[cl[j]] = rho_aver
                for i in range(ND):
                    if rho[i] < bord_rho[cl[i]]:
                        halo[i] = -1

            for i in range(NCLUST):
                nc = 0
                nh = 0
                for j in range(ND):
                    if cl[j] == i:
                        nc = nc + 1
                    if halo[j] == i:
                        nh = nh + 1
                print('CLUSTER: %i CENTER: %i ELEMENTS: %i CORE: %i HALO: %i' % (i + 1, icl[i] + 1, nc, nh, nc - nh))
            # print , start from 1

            ## save CLUSTER_ASSIGNATION
            print('Generated file:CLUSTER_ASSIGNATION')
            print('column 1:element id')
            print('column 2:cluster assignation without halo control')
            print('column 3:cluster assignation with halo control')
            clusters = np.array([np.arange(ND) + 1, cl + 1, halo + 1]).T
            np.savetxt('CLUSTER_ASSIGNATION_%.2f_%.2f_.txt' % (rhomin, deltamin), clusters, fmt='%d\t%d\t%d')
            print('Result are saved in file CLUSTER_ASSIGNATION_%.2f_%.2f_.txt' % (rhomin, deltamin))
            print('\n\nDrag the mouse pointer at a cutoff position in figure DECISION GRAPH and press   OR   Press key n to quit')
            ################# plot the data points with cluster labels
            cmap = cm.rainbow(np.linspace(0, 1, NCLUST))
            plot2(axarr, rho, delta, cmap, cl, icl, XY, NCLUST)
            plt.show()
            return ()

    while 1:
        print('\n\n在图中选择一个切断的点   OR   Press key n to quit')
        cid = f.canvas.mpl_connect('button_press_event', onclick)
        plt.show()  # get current figure
        nID = input()
        if nID == 'n':
            f.canvas.mpl_disconnect(cid)
            print('Saving the figure in file CLUSTER_ASSIGNATION.png')
            figure = plt.gcf()
            figure.set_size_inches(24, 8)
            figure.savefig('CLUSTER_ASSIGNATION.png', dpi=300)
            break


def mds(d, dimensions=2):
    """
    Multidimensional Scaling - Given a matrix of interpoint distances,
    find a set of low dimensional points that have similar interpoint
    distances.
    """
    print('Multidimensional Scaling for scatter plot...')
    (n, n) = d.shape
    E = (-0.5 * d ** 2)
    # Use mat to get column and row means to act as column and row means.
    Er = np.mat(np.mean(E, 1))
    Es = np.mat(np.mean(E, 0))
    F = np.array(E - np.transpose(Er) - Es + np.mean(E))
    [U, S, V] = svd(F)
    Y = U * np.sqrt(S)
    return (Y[:, 0:dimensions], S)


def readfile(file, dimensions=2, sep=' '):
    '''
    Input file format: 3 columns ,seperated by ' '
    Column 1: element1
    Column 2: element2
    Column 3: distance between element1 and element2
    For example: (id > 0)
        1   2   0.6
        1   3   2.3
        2   3   1.4
    Return (dist,xxdist,ND,N)
    '''
    print('Loading the file ...')
    xx = np.genfromtxt(file, delimiter=sep, names=['x', 'y', 'dist'], dtype="i8,i8,f8")
    # ND: number of data point
    X = xx['x']
    Y = xx['y']
    xxdist = xx['dist']
    ND = Y.max()
    NL = X.max()
    if NL > ND:
        ND = NL
    # N: number of point pairs/distance
    N = xx.shape[0]
    dist = np.zeros((ND, ND))

    # dist may save half of memory
    for i in range(N):
        ii = X[i] - 1
        jj = Y[i] - 1
        dist[ii, jj] = xxdist[i]
        dist[jj, ii] = xxdist[i]

    return ((dist, xxdist, ND, N))


def rhodelta(dist, xxdist, ND, N, percent=2.0):
    '''
    Input file format: 3 columns ,seperated by ' '
    Return (rho,delta,ordrho)
    '''
    print('Caculating rho and delta...')
    print('average percentage of neighbours (hard coded): %5.6f' % (percent))

    position = round(N * percent / 100)
    sda = np.sort(xxdist)
    dc = sda[position]
    print('Computing Rho with gaussian kernel of radius: %12.6f\n' % (dc))

    rho = np.zeros(ND)
    # Gaussian kernel
    for i in range(ND - 1):
        for j in range((i + 1), ND):
            rho[i] = rho[i] + np.exp(-(dist[i, j] / dc) * (dist[i, j] / dc))
            rho[j] = rho[j] + np.exp(-(dist[i, j] / dc) * (dist[i, j] / dc))

    maxd = dist.max()
    ordrho = (-rho).argsort()
    delta = np.zeros(ND)
    nneigh = np.zeros(ND)
    delta[ordrho[0]] = -1
    nneigh[ordrho[0]] = 0

    for ii in range(1, ND):
        delta[ordrho[ii]] = maxd
        for jj in range(ii):
            if dist[ordrho[ii], ordrho[jj]] < delta[ordrho[ii]]:
                delta[ordrho[ii]] = dist[ordrho[ii], ordrho[jj]]
                nneigh[ordrho[ii]] = ordrho[jj]

    nneigh1 = nneigh.astype(int)
    nneigh = nneigh1
    ordrho1 = ordrho.astype(int)
    ordrho = ordrho1

    delta[ordrho[0]] = delta.max()
    print('Generated file:DECISION GRAPH')
    print('column 1:Density')
    print('column 2:Delta\n')
    dg = np.array([rho, delta]).T
    np.savetxt('DECISION_GRAPH.txt', dg, fmt='%.4f')
    return ((rho, delta, ordrho, dc, nneigh))


def run(*args, **kwargs):
    '''
    return cluster id
    '''
    file = kwargs.get('fi')
    sep = kwargs.get('sep', ' ')
    ########
    (dist, xxdist, ND, N) = readfile(file, dimensions=2, sep=sep)
    XY, eigs = mds(dist)
    (rho, delta, ordrho, dc, nneigh) = rhodelta(dist, xxdist, ND, N, percent=2.0)
    DCplot(dist, XY, ND, rho, delta, ordrho, dc, nneigh, 17, 0.1)


filein = "test.dat"
run(fi=filein, sep='\t')
