import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

data = np.fromfile("test", dtype=np.int32).reshape((4096, 4096))
data = data.T
plt.figure(figsize=(10, 10))
plt.imshow(data, cmap="jet", norm=LogNorm())
plt.colorbar(label="Iteration Count")
plt.title("Mandelbrot Set (Log Scale)")
plt.show()