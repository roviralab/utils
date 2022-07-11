import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib import colors

rad2grad = 180/np.pi

def mercator(x,y,z,min_energia,max_energia):
    range_energy = max_energia
    #design specifications
    fig, ax = plt.subplots(figsize=(22,15))
    plt.rcParams.update({'font.size': 30})
    for tick in ax.get_xticklabels():
       tick.set_fontsize(28)
    for tick in ax.get_yticklabels():
       tick.set_fontsize(28)
    ax.axis((0, 359, 180, 0))
    plt.xticks(np.arange(0, 361, step=30))
    plt.yticks(np.arange(0, 181, step=45))
    ax.set_xlabel("φ (degrees)", fontsize=28)
    ax.set_ylabel("θ (degrees)", fontsize=28)

    #contorn lines for each energy kcal/mol
    levels = np.arange(0, range_energy, 1)
    ax.grid(color="k")

    #conformation labels

    phi_values = [0, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330, 360, 0, 30, 60, 90, 120, 150, 180,
                  210, 240, 270, 300, 330, 360, 0, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330, 360]
    theta_values = [54.7, 50.8, 54.7, 50.8, 54.7, 50.8, 54.7, 50.8, 54.7, 50.8, 54.7, 50.8, 54.7, 90, 90, 90, 90, 90,
                    90, 90, 90, 90, 90, 90, 90, 90,125.3, 129.2, 125.3, 129.2, 125.3, 129.2, 125.3, 129.2, 125.3, 129.2,
                    125.3, 129.2, 125.3]
    labels_conf = ["$^OE$", "$^OH_1$", "$E_1$", "$^2H_1$", "$^2E$", "$^2H_3$", "$E_3$", "$^4H_3$", "$^4E$", "$^4H_5$",
                   "$E_5$", "$^OH_5$", "$^OE$", "$^{3,O}B$", "$^3S_1$", "$B_{1,4}$", "$^5S_1$", "$^{2,5}B$", "$^2S_O$",
                   "$B_{3,O}$", "$^1S_3$","$^{1,4}B$", "$^1S_5$", "$B_{2,5}$", "$^OS_2$", "$^{3,O}B$", "$^3E$",
                   "$^3H_4$", "$E_4$", "$^5H_4$", "$^5E$", "$^5H_O$", "$E_O$", "$^1H_O$", "$^1E$", "$^1H_2$","$E_2$",
                   "$^3H_2$", "$^3E$"]
    ax.scatter(phi_values, theta_values, marker=".", color='black', s=400)
    for i, txt in enumerate(labels_conf):
        plt.annotate(txt, (phi_values[i], theta_values[i]-2), fontsize=38)
    # arrows and conformations names of the poles
    ax.annotate('', xy=(0,-0.09), xycoords='axes fraction', xytext=(1,-0.09), arrowprops=dict(arrowstyle='<->'))
    ax.annotate('', xy=(0,1.025), xycoords='axes fraction', xytext=(1,1.025), arrowprops=dict(arrowstyle='<->'))
    ax.text(0.42, 0.91, "$^4$C$_1$", transform=plt.gcf().transFigure, fontsize=38)
    ax.text(0.48, 0.005, "$^1$C$_4$", transform=plt.gcf().transFigure, fontsize=38)

    #matrix of energies for matplotlib to understand
    Z = z.reshape(int(np.sqrt(len(x))), int(np.sqrt(len(y))))
    T = Z.T #--> per quan cv1=phi i cv2=theta

    #colormap style
    bounds = np.linspace(0, range_energy, 1000, endpoint=True)
    cmap = cm.get_cmap('RdYlBu_r', 1000)
    norm = colors.BoundaryNorm(bounds, cmap.N)
    cmap = ax.imshow(T, aspect='auto', origin='lower', extent=(x.min(), x.max(), y.min(), y.max()),
                     cmap=cmap,
                     norm=norm
                     )

    #contorlines
    ax.contour(T, levels, colors='k', extent=(x.min(), x.max(), y.min(), y.max()))
    #legend colormap
    cbar = fig.colorbar(cmap, ax=ax, ticks=np.arange(0,range_energy,5), pad=0.07, aspect=16)
    cbar.set_label("Free Energy (kcal/mol)", fontsize=28, rotation=270, labelpad=34)

    fig.savefig("fes-mercator.png")

##MAIN##
input_file = "fes-mercator.dat"
pregunta = input("Is it reweighting (1) or sum_hills (2)? ")
if pregunta == "1":
    x, y, z, kk1, kk2 = np.loadtxt(input_file, unpack=True)
    x, y = x*rad2grad, y*rad2grad
    max_energia=30
    #substitute inf for the max value
    z = np.where(z == np.inf,max_energia,z)
    min_energia = min(z)
    mercator(y,x,z,min_energia,max_energia)
elif pregunta == "2":
    x, y, z = np.loadtxt(input_file, unpack=True)
    x, y = x * rad2grad, y * rad2grad
    min_energia = min(z)
    print(min_energia)
    z = z - min_energia
    max_energia = max(z)
    print(max_energia)
    mercator(x, y, z, min_energia, max_energia)

