import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib import colors

# rad2grad = 180/np.pi

def mercator(x,y):
    #design specifications
    fig, ax = plt.subplots(figsize=(20,15))
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

    ax.scatter(x,y, marker=".", c=range(len(x)),cmap="viridis", s=8000)

    fig.savefig("mercator_xtal_neg_noG2SGWT.png")

##MAIN##
input_file = "getpuck_all.out"
x,y=np.loadtxt(input_file, unpack=True)
mercator(y,x)
