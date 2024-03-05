import numpy as np
from collections import defaultdict
import matplotlib.pyplot as plt
from matplotlib import cm

import matplotlib.animation as anim

class population:
  def __init__(self, N, x, p):
    self.N = N # population density
    self.x = x # phenotype vector
    self.p = p # predation coef
    self.r = 1 - p**(0.9) # consumption coef
    self.d = len(x) # phenotype vector dimensions
  def __str__(self):
    return f"population density: {self.N}, phenotype vector: {self.x}, predation coef: {self.p}, phenotype vector dimensions: {self.d}"

  def mutate(self):
    N = 0.00001 # mutant population
    mu = 0.01   # standart mutation deviation
    x = self.x + np.random.normal(0, mu, size=2)
    p = self.p + np.random.normal(0, mu, size=1)
    p = p[0]
    if p<0:
      p=0
    if p>1:
      p=1
    new_p = population(N, x, p)
    return new_p

# CLUSTERING MODULE

# distance between two vectors
def dist(x, y, d):
  res = 0
  for i in range(d):
    res += (x[i]-y[i])**2
  res = np.sqrt(res)
  return res

# get one weighted population from the list of populations
def weight(N, poplist):
  coef = (poplist[0].N / N)
  x = poplist[0].x * coef
  p = poplist[0].p * coef
  for i in range(1, len(poplist)):
    coef = (poplist[i].N / N)
    x += poplist[i].x * coef
    p += poplist[i].p * coef
  
  new_pop = population(N, x, p)
  return new_pop

def pop_clust(poplist, eps):    # clustering population function:
                                # poplist - list of clusters, eps - min distance
  p = population(0, np.zeros(2), 0.0)
  vecp = [p]
  clust = defaultdict(lambda:vecp)
  i = len(poplist)-1            # go from the end
  while poplist[i].N < 1e-6/2:
    poplist.pop()
    i -= 1

  clust[0] = [poplist[i]]       # init first cluster
  L = len(poplist)
  cflag = False                 # if population does not contain in any cluster set True and add new cluster

  min_dist = 0.1
  key_closest = 0
  for i in range(L-2, -1, -1):
    if  poplist[i].N < 1e-6/2:   # if population is too smOll drop it
      poplist.pop()
    else:
      for key in clust.keys():
        d = dist(clust[key][0].x, poplist[i].x, poplist[i].d)
        if d < eps:
          if d<min_dist:
              min_dist = d
              key_closest = key
              cflag = True
              
      if cflag:
        clust[key_closest].append(poplist[i])     
      else:
        clust[len(clust.keys())] = [poplist[i]]
    cflag = False
    min_dist = 0.1
    key_closest = 0

  new_poplist = []
  for key in clust.keys():
    N = clust[key][0].N
    for i in range(1, len(clust[key])):
      N += clust[key][i].N
    p = weight(N, clust[key]) 
    new_poplist.append(p)

  return new_poplist

def pop_clust2(poplist, eps): 
  i = len(poplist)-1            # go from the end
  while poplist[i].N < 1e-6/2:
    poplist.pop()
    i -= 1
  newp = [poplist[0]]
  if (len(poplist)>1):
      min_dist = 0.1
      cflag = False
      l = len(poplist)
      bufp = poplist
      flag = np.full(l, True)
      newp = []
      for j in range(l):
          if flag[j]:
              for i in range(l):
                  if flag[i]:
                      d = dist(bufp[i].x, bufp[j].x, 2)
                      if (d<eps and abs(bufp[j].p-bufp[i].p)<0.1 and i != j):
                          if (d<min_dist):
                              min_dist=d
                              key_closest = i
                          cflag=True
          if cflag:
              w = bufp[j].N + bufp[key_closest].N
              bufp[key_closest].p = (bufp[key_closest].p*bufp[key_closest].N + bufp[j].p*bufp[j].N)/w
              bufp[key_closest].x = (bufp[key_closest].x*bufp[key_closest].N + bufp[j].x*bufp[j].N)/w
              bufp[key_closest].N = w
              flag[j] = False
          
          cflag=False
          min_dist = 0.1
      for j in range(l):
         if flag[j]:
             newp.append(bufp[j])
  return newp
    
    
# Output video writer
FFMpegWriter = anim.writers['ffmpeg']
metadata = dict(title='Movie Test', artist='Matplotlib', comment='Movie support!')
writer = FFMpegWriter(fps=64, metadata=metadata)
fig = plt.figure(figsize=(6.4, 4.8), dpi=200)


inpass = ''
outpass = ''
ch = 0
if ch==0:
    inpass = "diploma_refactored/populationz"
    outpass = "diploma_refactored/video.mp4"
    
    outpass2 = "diploma_refactored/qdynamic.png"
elif ch==1:
    inpass = "diploma_refactored/comp1/populationz"
    outpass = "diploma_refactored/video_comp1.mp4"
    
    outpass2 = "diploma_refactored/qdynamic_comp1.png"
else:
    inpass = "diploma_refactored/comp2/populationz"
    outpass = "diploma_refactored/video_comp2.mp4"
    
    outpass2 = "diploma_refactored/qdynamic_comp2.png"

writer.setup(fig, outfile=outpass, dpi=200)

for i in range(0, 2500501, 500):
    # Open file with populations
    f = open(f"{inpass}/{i}.txt", 'r')
    #it = int(f.name.split('.')[0])
    lines = f.read()
    lines = lines.split('\n')
    poplist = []
    for line in lines:
      line = line.split()
      p = population(float(line[0]), np.array(line[1:-1]).astype(float), float(line[-1]))
      poplist.append(p)
    f.close()
    
    view_poplist = pop_clust2(poplist, 0.08)

    # Визуализация точек
    fig.clear()
    ax = plt.gca()
    cmap = cm.get_cmap('Spectral_r')  # Создаем цветовую карту
    norm = plt.Normalize(0, 1)     # Нормализуем значения для цветовой карты

    for point in view_poplist:
        circle1 = plt.Circle((point.x[0], point.x[1]), np.sqrt(point.N)/30.0, color=cmap(norm(point.p)), alpha=0.3)
        plt.gca().add_patch(circle1)

    # Добавляем colorbar
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    plt.colorbar(sm)
    ax.set_axisbelow(True)
    ax.set_xlim([-2, 2])
    ax.set_ylim([-2, 2])
    plt.xlabel("phenotype 1",fontsize=7)
    plt.ylabel("phenotype 2",fontsize=7)
    plt.grid(True, linestyle='--')
    plt.title(f"Population dynamic. Time unit = {i}")
    #plt.savefig(f"{it}.png", dpi=300, bbox_inches='tight')
    writer.grab_frame()

    if i%100000==0:
        print(f"{i} iterations have passed")

writer.finish()