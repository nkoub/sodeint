import numpy as np
import matplotlib.pyplot as plt
sp = np.loadtxt("ring_N100_D0.0000_coupling0.08000_states.out")
plt.imshow(sp,interpolation="nearest",origin="upper",cmap=plt.cm.binary,aspect="auto")
plt.show()

