
import numpy as np
from sklearn import datasets
import matplotlib.pyplot as plt
import matplotlib as mpl
from neupy import algorithms, utils
mintemp=33.2; maxtemp=34.85; division=0.05
cmap = mpl.cm.jet #gist_rainbow_r
norm = mpl.colors.Normalize(vmin=mintemp, vmax=maxtemp)

plt.style.use('ggplot')
utils.reproducible()

k=0
classify = []
for line in open('Array_temp.txt','r').readlines()[:]:
	classify.append([float(i) for i in line.split( )])
	k = k+1
	print(k)
"""
	a = [float(i) for i in line.split( )]
	fr = int(len(a)/3)
	a1 = a[0:fr]
	a2 = a[fr:fr*2]
	a3 = a[fr*2:]
	classify.append(list([a1,a2,a3]))
	k = k+1
	print(k)
	temp.append(float(line.split( )[0]))
	salt.append(float(line.split( )[1]))
	zeta.append(float(line.split( )[2]))
	mag.append(float(line.split( )[3]))
	ang.append(float(line.split( )[4]))
	classify.append(list([salt[-1],mag[-1],ang[-1],zeta[-1]]))
"""
classify = np.array(classify)

ggplot_colors = plt.rcParams['axes.prop_cycle']
colors = np.array([c['color'] for c in ggplot_colors])

#dataset = datasets.load_iris()
# use only two features in order
# to make visualization simpler
# data = dataset.data[:, [2, 3]]
data = classify[:]
#target = kk #dataset.target

#dataset = datasets.load_iris()
#data = dataset.data[:, [2, 3]]
#target = dataset.target

sofm = algorithms.SOFM(
# Use only two features for the input
n_inputs=7320,

# Number of outputs defines number of features
# in the SOFM or in terms of clustering - number
# of clusters
n_outputs=4,

# In clustering application we will prefer that
# clusters will be updated independently from each
# other. For this reason we set up learning radius
# equal to zero
learning_radius=0,

# Training step size or learning rate
step=0.25,

# Shuffles dataset before every training epoch.
shuffle_data=True,

# Instead of generating random weights
# (features / cluster centers) SOFM will sample
# them from the data. Which means that after
# initialization step 3 random data samples will
# become cluster centers
weight='sample_from_data',

# Shows training progress in terminal
verbose=True,
)
sofm.train(data[0:15*12], epochs=200)

n_outputs = 4
fil = open('Nodes_temp.txt', 'w')
for i in range(n_outputs):
	arr = [str(i) for i in np.round(sofm.weight[:,i],1)]
	fil.write('  '.join(arr) + '\n')
fil.close()

p1 = sofm.predict(data[0:15*12])
p2 = sofm.predict(data[15*12:])

arr1,arr2 = [],[]
for i in range(len(p1)):
	arr1.append(np.where(p1[i]==1)[0][0])
for i in range(len(p2)):
	arr2.append(np.where(p2[i]==1)[0][0])
#print('Case 1 :', arr1.count(0)/(len(arr1)+1),'Case 2 :',arr1.count(1)/(len(arr1)+1),'Case 3 :',arr1.count(2)/(len(arr1)+1))
#print('Case 1 :', arr2.count(0)/(len(arr2)+1),'Case 2 :',arr2.count(1)/(len(arr2)+1),'Case 3 :',arr2.count(2)/(len(arr2)+1))

month1,month2 = [[] for i in range(12)],[[] for i in range(12)]
for i in range(len(arr1)):
	month1[int(i%12)].append(arr1[i])
for i in range(len(arr2)):
	month2[int(i%12)].append(arr2[i])

m1,m2 = [[],[],[],[]],[[],[],[],[]]
for i in range(12):
	m1[0].append(100*month1[i].count(0)/len(month1[1]))
	m1[1].append(100*month1[i].count(1)/len(month1[1]))
	m1[2].append(100*month1[i].count(2)/len(month1[1]))
	m1[3].append(100*month1[i].count(3)/len(month1[1]))

	m2[0].append(100*month2[i].count(0)/len(month2[1]))
	m2[1].append(100*month2[i].count(1)/len(month2[1]))
	m2[2].append(100*month2[i].count(2)/len(month2[1]))
	m2[3].append(100*month2[i].count(3)/len(month2[1]))

m1,m2 = np.array(m1),np.array(m2)

colors = ['b','m','y','c']
fig, ax = plt.subplots()
for i in range(len(arr1)):
	plt.scatter(np.mean(data[i,3720:]),np.mean(data[i,0:3660]), color = colors[int(arr1[i])])

for j in range(4):
	plt.scatter(np.mean(sofm.weight[3660:,j]),np.mean(sofm.weight[0:3720,j]), s=250, color = colors[j], edgecolor='black')
ax.set_xlabel('Southern Hemisphere [C]')
ax.set_ylabel('Norhtern Hemisphere [C]')
fig.savefig("./Cluster_map.png" , dpi=200, bbox_inches='tight')

fig, ax = plt.subplots()
ind = np.arange(12)
width = 1
p1 = plt.bar(ind, m1[0], width, color = 'b', label='Winter')
p2 = plt.bar(ind, m1[1], width,bottom=m1[0], color = 'm', label='Transition (Warm)')
p2 = plt.bar(ind, m1[2], width,bottom=m1[0]+m1[1], color = 'y', label='Summer')
p2 = plt.bar(ind, m1[3], width,bottom=m1[0]+m1[1]+m1[2], color = 'c', label='Transition (Cold)')
ax.set_xticks(ind + width / 2 - 0.5)
ax.set_xticklabels(('Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'))
plt.xlim(-0.5,11.5)
plt.ylim(0,100)
ax.set_ylabel('Frequency')


# Shrink current axis's height by 10% on the bottom
box = ax.get_position()
ax.set_position([box.x0, box.y0 + box.height * 0.1,
                 box.width, box.height * 0.9])

# Put a legend below current axis
ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.1),
          fancybox=True, shadow=True, ncol=4)


fig.savefig("./Seasons_1979_1991.png" , dpi=200, bbox_inches='tight')


fig, ax = plt.subplots()
ind = np.arange(12)
width = 1
p1 = plt.bar(ind, m2[0], width, color = 'b', label='Winter')
p2 = plt.bar(ind, m2[1], width,bottom=m2[0], color = 'm', label='Transition (Warm)')
p2 = plt.bar(ind, m2[2], width,bottom=m2[0]+m2[1], color = 'y', label='Summer')
p2 = plt.bar(ind, m2[3], width,bottom=m2[0]+m2[1]+m2[2], color = 'c', label='Transition (Cold)')
ax.set_xticks(ind + width / 2 - 0.5)
ax.set_xticklabels(('Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'))
plt.xlim(-0.5,11.5)
plt.ylim(0,100)
ax.set_ylabel('Frequency')


# Shrink current axis's height by 10% on the bottom
box = ax.get_position()
ax.set_position([box.x0, box.y0 + box.height * 0.1,
                 box.width, box.height * 0.9])

# Put a legend below current axis
ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.1),
          fancybox=True, shadow=True, ncol=4)


fig.savefig("./Seasons_1991_2017.png" , dpi=200, bbox_inches='tight')

