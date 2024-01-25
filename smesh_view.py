import numpy as np
import matplotlib.pyplot as plt
from matplotlib import patches
import matplotlib as mpl
import matplotlib.colors as colors
from mpl_toolkits.axes_grid1 import make_axes_locatable

mpl.rcParams.update({"font.size": 8})

def load_fortran_binary(filename):
    # Fortran data structure:
    #     4 byte head --> N
    #     N byte first record (int32, the shape of the object contained in the second record)
    #     4 byte foot == head
    #     4 byte head --> M
    #     M byte second record (int32 or float64, the matrix or vector)
    #     4 byte foot == head
    with open(filename) as f:
        # read shape
        record_length = np.fromfile(f, dtype=np.int32, count=1, sep="")[0]
        shape = np.fromfile(f, dtype=np.int32, count=record_length//4, sep="")
        dummy = np.fromfile(f, dtype=np.int32, count=1, sep="")
        # compute matrix size
        nelem = np.prod(shape)
        # here we use the record length to determine the type of data (8bit/element -> real, 4bit/element -> integer)
        record_length = np.fromfile(f, dtype=np.int32, count=1, sep="")[0]
        # select data type
        if record_length//nelem == 4:
            dtype = np.int32
        else:
            dtype = np.float64
        # read data
        if shape.size > 1:
            return np.fromfile(f, dtype=dtype, count=nelem, sep="").reshape(shape[::-1]).transpose()
        else:
            return np.fromfile(f, dtype=dtype, count=nelem, sep="").flatten()

class Triangulation():

    def __init__(self, pt, ve, vor_pt, vor_ve, vor_veb):
        self.pt = pt
        self.ve = ve
        self.vpt = vor_pt
        self.vve = vor_ve
        self.vveb = vor_veb
        self.kp = self.ve.shape[0]
        self.kpv = self.vveb.shape[0]
        index = np.zeros(pt.shape[0])
        for i in range(self.kp):
            index[self.ve[i,:]] = index[self.ve[i,:]] + i
        self.index = index/3.0
        nelem = ve.shape[0]

    def display(self, boundaries=True, plot_index=False, show_delaunay=True, show_voronoi=False, individual_voronoi=False):
        fig = plt.figure()
        ax = [fig.add_subplot(1,1,1)]
        if plot_index:
            ax[0].tripcolor(self.pt[:,0], self.pt[:,1], self.index, linewidths=0.0, cmap=plt.get_cmap("viridis"))
        nan_nan = np.array([np.nan, np.nan])
        if show_delaunay:
            print("building segments")
            segments = np.vstack([
                np.vstack([self.pt[self.ve[i,:]], self.pt[self.ve[i,0]], nan_nan]) for i in range(self.kp)
                ])
        if show_voronoi:
            print("building segments")
            segments_voronoi = np.vstack([
                np.vstack([self.vpt[self.vve[self.vveb[i,0]:self.vveb[i,1]]], self.vpt[self.vve[self.vveb[i,0]]], nan_nan]) for i in range(self.kpv)
                ])
        # plot all segments at once
        color = "black" if plot_index else "darkorange"
        alpha = 0.5 if plot_index else 1
        if show_delaunay:
            print("plotting segments")
            ax[0].plot(segments[:,0], segments[:,1], color=color, alpha=alpha, lw=0.2)
        if show_voronoi:
            print("plotting segments")
            if individual_voronoi:
                for i in range(self.kpv):
                    a = np.vstack([self.ptv[self.vev[i,:self.nvev[i]]], self.ptv[self.vev[i,0]]])
                    x = a[:,0]
                    y = a[:,1]
                    cx = np.mean(x)
                    cy = np.mean(y)
                    omega = 0.95
                    x = omega*x + (1 - omega)*cx
                    y = omega*y + (1 - omega)*cy
                    ax[0].plot(x, y, color="dodgerblue")
            else:
                ax[0].plot(segments_voronoi[:,0], segments_voronoi[:,1], color="dodgerblue", alpha=alpha, lw=0.2)
        xmax, xmin = np.amax(self.pt[:,0]), np.amin(self.pt[:,0])
        ymax, ymin = np.amax(self.pt[:,1]), np.amin(self.pt[:,1])
        bx = xmax - xmin
        by = ymax - ymin
        margin = 0.05*min(bx, by)
        for a in ax:
            a.set_xlim(xmin - margin, xmax + margin)
            a.set_ylim(ymin - margin, ymax + margin)
            a.set_aspect("equal")
        plt.tight_layout(pad=0.5, w_pad=0, h_pad=0)
        # plt.get_current_fig_manager().window.state("zoomed")
        plt.show()

def fast_figure(x, y):
    fig = plt.figure()
    ax = [fig.add_subplot(1, 1, 1)]
    xmax, xmin = np.amax(x), np.amin(x)
    ymax, ymin = np.amax(y), np.amin(y)
    bx = xmax - xmin
    by = ymax - ymin
    margin = 0.05*min(bx, by)
    # ax[0].set_xlim(xmin - margin, xmax + margin)
    # ax[0].set_ylim(ymin - margin, ymax + margin)
    ax[0].set_aspect("equal")
    plt.tight_layout(pad=0.5, w_pad=0, h_pad=0)
    # plt.get_current_fig_manager().window.state("zoomed")
    return fig, ax

fast = False

if fast:
    pt = load_fortran_binary(r"output_data\dt_pt.dat")
    ed = load_fortran_binary(r"output_data\dt_edge.dat")
    ed = ed - 1
    x = pt[0,:]
    y = pt[1,:]
    sx = np.column_stack([x[ed[0,:]], x[ed[1,:]], np.nan*np.zeros(ed.shape[1])]).flatten()
    sy = np.column_stack([y[ed[0,:]], y[ed[1,:]], np.nan*np.zeros(ed.shape[1])]).flatten()
    fig, ax = fast_figure(x, y)
    ax[0].plot(sx, sy, color="darkorange", alpha=1, lw=0.3) 
    plt.show()

else:
    pt          = load_fortran_binary(r"output_data\dt_pt.dat")
    ve          = load_fortran_binary(r"output_data\dt_ve.dat")
    vor_pt      = load_fortran_binary(r"output_data\dt_vor_pt.dat")
    vor_ve      = load_fortran_binary(r"output_data\dt_vor_ve.dat")
    vor_veb     = load_fortran_binary(r"output_data\dt_vor_veb.dat")
    ve = ve - 1
    vor_ve = vor_ve - 1
    vor_veb[0,:] = vor_veb[0,:] - 1
    pt          = np.transpose(pt)
    ve          = np.transpose(ve)
    vor_pt      = np.transpose(vor_pt)
    vor_veb     = np.transpose(vor_veb)

    print("Triangulation loaded.")

    triangulation = Triangulation(pt, ve, vor_pt, vor_ve, vor_veb)
    triangulation.display(show_delaunay=True, show_voronoi=True)
