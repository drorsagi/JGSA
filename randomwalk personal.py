

def calcnoise(lst): # calculates the noise for each metabolite
    ans = statistics.pstdev(lst)/statistics.mean(lst)
    return ans



import statistics
import numpy as np
import matplotlib.pyplot as plt
from numpy.random import randint, normal, uniform

# simulating our case: 400 metabolites, 500 steps, 20 walkers, each with individual step
# need to chnage here - each individual needs different pace

# # making all three time points on the same graph - 100, 200, 500
timepoints0 = [100, 150, 200]
n0 = 20 # animas
m0 = 434 # features
noiselst0 = [] # each element is a list of al noises for a particular time point
for ii in timepoints0:
    t_max = ii

    valuesanimal = []
    for jj in range(0, n0):  # for each animal we generate the metabolomics table
        tmp = (randint(0, 20)) / 100  # allows 20% variance per step
        # print(tmp) looking good
        random_numbers = randint(0, 2, size=(t_max, m0))  # a table of 500 steps and 434 metabolite per animal, each met flactuates with random walk
        steps = np.where(random_numbers == 0, -1 - tmp, +1 + tmp)
        # print(steps[0, 0]) - yey, adds the individual tmp to +- 1
        values0 = np.sum(steps, axis=0)  # for each metabolite we have the distance
        values1 = [ii + 1000 for ii in values0]
        valuesanimal.append(values1)  # each element is an animal
    print(len(valuesanimal), len(valuesanimal[0]))  # 20 434

    # calculating noise
    valuesanimalarr = np.array(valuesanimal)
    valuesfornoise = valuesanimalarr.transpose()
    print(len(valuesfornoise), len(valuesfornoise[0]))  # 434 20
    noise0 = []
    for kk in valuesfornoise:
        tmp = calcnoise(kk)
        noise0.append(tmp)
    print(len(noise0)) # 434
    noiselst0.append(noise0)


# plotting the histogram
medyoung = statistics.median(noiselst0[0])
medold = statistics.median(noiselst0[1])
medal = statistics.median(noiselst0[2])

fig, ax = plt.subplots()
bins = np.linspace(0, 0.03, 50)

plt.hist(noiselst0[0], bins, alpha=0.5, color = 'blue', label='Young')
plt.hist(noiselst0[1], bins, alpha=0.5, color = 'red', label='old')
plt.hist(noiselst0[2], bins, alpha=0.5, color = 'green', label='AL')

plt.axvline(x=medyoung, ymin=0, ymax=1, color = 'blue', linewidth = 4)
plt.axvline(x=medold, ymin=0, ymax=1, color = 'red', linewidth = 4)
plt.axvline(x=medal, ymin=0, ymax=1, color = 'green', linewidth = 4)

ax.tick_params(axis='x', labelsize=22)
ax.tick_params(axis='y', labelsize=22)
# fig.savefig('random walk personal.pdf')
plt.show()

# plotting the path of 5 individuals - they become wider with time -correct for individual pace!!!
t_max = 100
n0 = 5 # 10
animalval0 = [] #
for ii in range(0, n0):
    tmp = (randint(0, 150)) / 100  # allows 150% variance per step
    random_numbers = randint(0, 2, size=t_max)
    steps = np.where(random_numbers == 0, -1 - tmp, +1 + tmp)
    values = np.cumsum(steps)
    valueslst = values.tolist()
    animalval0.append([0] + valueslst) # initial value is set to 0
print(len(animalval0), len(animalval0[0]), animalval0[0][0], animalval0[0][1]) # 5 101 0 1.08 - yey

colors0 = ['red', 'blue', 'green', 'black', 'pink']
markers0 = ['o', 'v', 'p', '*', '<']
t00 = [ii for ii in range(0, t_max+1)]
fig, ax = plt.subplots()
for ii in range(0, n0):
    plt.plot(t00, animalval0[ii], color = colors0[ii], marker = markers0[ii], linestyle='dashed') # for each row works, how to plot in paralel all rows?
ax.tick_params(axis='x', labelsize=22)
ax.tick_params(axis='y', labelsize=22)
# fig.savefig('random walk personal path.pdf')
plt.show() #






