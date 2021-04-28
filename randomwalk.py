
# https://sphelps.net/teaching/scf/slides/random-walks-slides.html



def calcnoise(lst): # calculates the noise for each metabolite, adjusting for mean = 0
    if statistics.mean(lst) != 0:
        ans = statistics.pstdev(lst)/statistics.mean(lst)
    elif statistics.mean(lst) == 0:
        ans = 0
    return ans


import statistics
import numpy as np
import matplotlib.pyplot as plt
from numpy.random import randint, normal, uniform



# simulating our case: 400 metabolites, 1000 steps, 20 walkers


# making all three time points on the same graph - 100, 200, 500
timepoints0 = [100, 150, 200]
n = 20 # animas
m0 = 434 # features
noiselst0 = [] # each element is a list of al noises for a particular time point
for ii in timepoints0:
    t_max = ii

    valuesmet0 = []
    for jj in range(0, m0):
        random_numbers = randint(0, 2, size=(t_max, n))
        steps = np.where(random_numbers == 0, -1, +1)
        values0 = np.sum(steps, axis=0) # for each animal we have the distance
        values1 = [ii + 1000 for ii in values0]
        valuesmet0.append(values1) # each element is a metabolite
    print(len(valuesmet0), len(valuesmet0[0])) # 434 20 - yey

    # calculating noise
    noise0 = []
    for kk in valuesmet0:
        tmp = calcnoise(kk)
        noise0.append(tmp)
    print(len(noise0))
    noiselst0.append(noise0)


medyoung = statistics.median(noiselst0[0]) # + 0.1
medold = statistics.median(noiselst0[1]) # + 0.1
medal = statistics.median(noiselst0[2]) # + 0.1

fig, ax = plt.subplots()
bins = np.linspace(0, 0.03, 75)

plt.hist(noiselst0[0], bins, alpha=0.5, color = 'blue', label='Young')
plt.hist(noiselst0[1], bins, alpha=0.5, color = 'red', label='old')
plt.hist(noiselst0[2], bins, alpha=0.5, color = 'green', label='AL')


plt.axvline(x=medyoung, ymin=0, ymax=1, color = 'blue', linewidth = 4)
plt.axvline(x=medold, ymin=0, ymax=1, color = 'red', linewidth = 4)
plt.axvline(x=medal, ymin=0, ymax=1, color = 'green', linewidth = 4)

ax.tick_params(axis='x', labelsize=22)
ax.tick_params(axis='y', labelsize=22)
# fig.savefig('random walk.pdf')
plt.show()

# plotting the path of 5 individuals - they become wider with time
t_max = 100
n = 5 # 10
random_numbers = randint(0, 2, size=(t_max, n))
steps = np.where(random_numbers == 0, -1, +1)
values = np.cumsum(steps, axis=0)
print(values)
valueslst = values.tolist()
initval0 = [[0]*n]
print(initval0)
values1 = initval0 + valueslst
valuesarrt = np.array(values1).transpose()
colors0 = ['red', 'blue', 'green', 'black', 'pink']
markers0 = ['o', 'v', 'p', '*', '<']
t00 = [ii for ii in range(0, t_max+1)]
fig, ax = plt.subplots()
for ii in range(0, n):
    plt.plot(t00, valuesarrt[ii], color = colors0[ii], marker = markers0[ii], linestyle='dashed') # for each row works, how to plot in paralel all rows?
ax.tick_params(axis='x', labelsize=22)
ax.tick_params(axis='y', labelsize=22)
plt.show() #











