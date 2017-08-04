import matplotlib
# fix the display settings for the cluster
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import numpy as np
import pyemma
from matplotlib.pyplot import *

Rg = np.genfromtxt('Rg.xvg')[:,1]
Qhel = np.genfromtxt('Qhel.dat')

H_Rg_Qhel, Rg_edges, Qhel_edges = np.histogram2d(Rg, Qhel, bins=50, range=None, normed=None, weights=None)

Rg_bins = 0.5*(Rg_edges[1:] + Rg_edges[:-1])
Qhel_bins = 0.5*(Qhel_edges[1:] + Qhel_edges[:-1])

fignum=0
fig = plt.figure(fignum)

#mycmap_FES = cm.get_cmap('gist_earth')
#mycmap_FES = cm.get_cmap('OrRd_r')
#mycmap_FES = cm.get_cmap('BuPu_r')
#mycmap = cm.get_cmap('ocean')
mycmap_FES = cm.get_cmap('terrain')
#mycmap.set_over('w')
#_jet_data = {'red': ((0., 0, 0), (0.35, 0, 0), (0.66, 1, 1), (0.89,1, 1), (1, 0.5, 0.5)), 'green': ((0., 0, 0), (0.125,0, 0), (0.375,1, 1), (0.64,1, 1), (0.91,0,0), (1, 0, 0)), 'blue': ((0., 0.5, 0.5), (0.11, 1, 1), (0.34, 1, 1), (0.65,0, 0), (1, 0, 0))}
#cdict = {'red': ((0., 0, 0), (0.35, 0, 0), (0.66, 1, 1), (0.89,1, 1), (1., 1, 1)), 'green': ((0., 0, 0), (0.125,0, 0), (0.375,1, 1), (0.64,1, 1), (0.91,0,0), (1., 1, 1)), 'blue': ((0., 0.5, 0.5), (0.11, 1, 1), (0.34, 1, 1), (0.65,0, 0), (1., 1, 1))}
#mycmap = matplotlib.colors.LinearSegmentedColormap('my_colormap',cdict,256)

H_Rg_Qhel = np.rot90(H_Rg_Qhel)
H_Rg_Qhel = np.flipud(H_Rg_Qhel)

FE_H = -np.ma.log(H_Rg_Qhel)
FE_H -= np.min(FE_H)
FE_max = np.max(FE_H)
#FE_H.set_fill_value(np.max(FE_H))
FE_H = FE_H.filled()
plt.pcolormesh(Rg_bins, Qhel_bins, FE_H, cmap='terrain',vmax=FE_max)
#plt.scatter(Rg, Qhel)


#cbar = plt.colorbar()
#cbar.set_ticks([0,3,6])
#cbar.ax.set_yticklabels([r'0', r'3', r'6'],fontweight='normal')
#cbar.ax.tick_params(labelsize='32',width=4.0,length=21)
#cbar.set_label(r'$k_{\rm B} T$', rotation=270, fontsize='40', labelpad=35)
#cbar.ax.get_children()[2].set_linewidth(4.0)
#cbar.ax.get_children()[3].set_linewidth(4.0)
#cbar.update_ticks()

#myColorMap = jet;
#myColorMap(1,:) = [1 1 1];
#colormap(myColorMap);
#colorbar;
#ma=max(value);
#mi=min(value);

#plt.axes().set_aspect(305)
ax = fig.add_subplot(111)
#plt.xlabel(r'$\Psi$ (deg)',fontsize='44',fontweight='normal')
#plt.ylabel(r'$R_{\operatorname{\boldsymbol{1-4}}}$ (nm)',fontsize='44',fontweight='normal',labelpad=7)
#plt.xticks([-180.0,-90.0,0,90.0,180.0],['-180','','0','','180'],fontsize='32',fontweight='normal')
#plt.yticks([0.4,0.6,0.8,1.0],['0.4','0.6','0.8','1.0'],fontsize='32',fontweight='normal')
#ax.tick_params(axis='both', which='major', pad=10)

#PAX = plt.gca()
#PAX.yaxis.set_label_coords(-0.04, 0.5)
#PAX.xaxis.set_label_coords(0.5, -0.06)
#plt.xlim([-180,180])
#plt.ylim([0.33,1.075])
#plt.legend()
#legend = plt.legend(fontsize='10',loc='upper center', bbox_to_anchor=(0.15, 1.5))
#frame = legend.get_frame()
#frame.set_linewidth('3.0')
#ltext = legend.get_texts()
#plt.setp(ltext, fontweight='bold')
ax.spines['left'].set_linewidth(5.0)
ax.spines['right'].set_linewidth(5.0)
ax.spines['top'].set_linewidth(5.0)
ax.spines['bottom'].set_linewidth(5.0)


plt.tight_layout()
plt.savefig('fig_FES_Rg_Qhel.png')
plt.close(fig)
fignum += 1

