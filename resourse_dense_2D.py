import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

import matplotlib.animation as anim

# Output video writer
FFMpegWriter = anim.writers['ffmpeg']
metadata = dict(title='Movie Test', artist='Matplotlib', comment='Movie support!')
writer = FFMpegWriter(fps=64, metadata=metadata)
fig = plt.figure(figsize=(6.4, 4.8), dpi=200)

inpass = ''
outpass = ''
ch = 0
if ch==0:
    inpass = "diploma_refactored/res_dir"
    outpass = "diploma_refactored/video_res_2D.mp4"
elif ch==1:
    inpass = "diploma_refactored/comp1/res_dir"
    outpass = "diploma_refactored/video_comp1_res_2D.mp4"
else:
    inpass = "diploma_refactored/comp2/res_dir"
    outpass = "diploma_refactored/video_comp2_res_2D.mp4"

writer.setup(fig, outfile=outpass, dpi=200)
x = np.linspace(-3, 3, 60)
y = np.linspace(-3, 3, 60)
xx, yy = np.meshgrid(x, y)
xs, ys = np.meshgrid(x, y, sparse=True)
R0 = np.full([xs.shape[1],ys.shape[0]], 1e-05)
for it in range(0, 2500501, 500):
    # Open file with populations
    with open(f"{inpass}/{it}.txt", "r") as myfile:
        info = myfile.read()
        info = info.strip()
        res_distribution = info.split(' ')
        for i in range(60):
            for j in range(60):
                R0[j][i] = float(res_distribution[i*60+j])

    fig.clear()
    ax = plt.gca()
    cmap = cm.get_cmap('inferno')  # Создаем цветовую карту
    norm = plt.Normalize(0, 50)     # Нормализуем значения для цветовой карты
    ax.set_xlim([-3, 3])
    ax.set_ylim([-3, 3])
    h = plt.contourf(x, y, R0, levels=100, cmap='inferno', norm=norm)
    plt.axis('scaled')
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    plt.colorbar(sm)
    plt.title(f"Resource density 2D. Time unit = {it}")
    writer.grab_frame()
    if it%100000==0:
        print(f"{it} iterations have passed")

writer.finish()