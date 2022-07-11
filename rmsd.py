import numpy as np
import matplotlib.pyplot as plt

input_file1 = "rmsd.dat"
inp=np.loadtxt(input_file1)

#figure
fig, ax = plt.subplots(figsize=(20,10))
for tick in ax.get_xticklabels():
    tick.set_fontsize(32)
for tick in ax.get_yticklabels():
    tick.set_fontsize(32)
ax.set_ylabel("rmsd ($\AA$)", fontsize=40)
ax.set_xlabel("Simulation time (ps)", fontsize=40)
ax.grid()
colors = ["#88E0EF", "#161E54", "#FF5151", "#FF9B6A","#88E0EF", "#161E54", "#FF5151", "#FF9B6A"]

ax.plot(inp[:,0]/100, inp[:,1], color=colors[0], linewidth=3.0, alpha=1, label="protein")
ax.plot(inp[:,0]/100, inp[:,2], color=colors[1], linewidth=3.0, alpha=0.75, label="substrate")
ax.plot(inp[:,0]/100, inp[:,3], color=colors[2], linewidth=3.0, alpha=0.75, label="active site")
ax.plot(inp[:,0]/100, inp[:,4], color=colors[3], linewidth=3.0, alpha=0.75, label="W434,W286")

fig.legend(fontsize=36, ncol=2)
plt.show()
fig.savefig("rmsd.png")
