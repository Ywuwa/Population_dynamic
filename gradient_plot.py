import numpy as np
from collections import defaultdict
import matplotlib.pyplot as plt



inpass = ''
outpass = ''
ch = 0
if ch==0:
    inpass = "diploma_refactored/populationz/selective_gradient_dynamic.txt"
    outpass = "diploma_refactored/pdynamic.png"
    outpass2 = "diploma_refactored/qdynamic.png"
    outpass3 = "diploma_refactored/gdynamic.png"
elif ch==1:
    inpass = "diploma_refactored/comp1/populationz/selective_gradient_dynamic.txt"
    outpass = "diploma_refactored/pdynamic_comp1.png"
    outpass2 = "diploma_refactored/qdynamic_comp1.png"
    outpass3 = "diploma_refactored/gdynamic_comp1.png"
else:
    inpass = "diploma_refactored/comp2/populationz/selective_gradient_dynamic.txt"
    outpass = "diploma_refactored/pdynamic_comp2.png"
    outpass2 = "diploma_refactored/qdynamic_comp2.png"
    outpass3 = "diploma_refactored/gdynamic_comp2.png"

f = open(inpass, 'r')
lines = f.read()
lines = lines.strip()
lines = lines.split(' ')
f.close()

gdynamic = [float(i) for i in lines]


x = np.arange(0, len(gdynamic)*500, 500)

#x = np.arange(1.0/len(gdynamic)*500, 1, 1.0/500)
plt.figure(figsize=(6.4, 4.8))
plt.gca().set(xlim=(0, max(x)), ylim=(0.0, max(gdynamic)))

#plt.margins(x=0, y=0)
plt.plot(x, gdynamic)
#plt.title('')
plt.xlabel('время')
plt.ylabel('сумма длин градиентов')

plt.savefig(outpass3, dpi=300, bbox_inches='tight')