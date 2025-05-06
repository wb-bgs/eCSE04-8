#!/usr/bin/env python
  
import numpy as np
import matplotlib.pyplot as plt


accuracy = 4


energies = {}
energies['cpu'] = {}
energies['gpu'] = {}

energies['cpu'][ '200'] = np.array([     44255,     44352,     44183 ])
energies['cpu'][ '300'] = np.array([    269873,    276444,    275851 ])
energies['cpu'][ '720'] = np.array([   7361281,   7364390,   7376685 ])
energies['cpu']['1440'] = np.array([ 140369359, 139858912, 141763753 ])

energies['gpu'][ '200'] = np.array([     14247,     15103,     14192 ])
energies['gpu'][ '300'] = np.array([     88466,     87348,     86453 ])
energies['gpu'][ '720'] = np.array([   3478627,   3435897,   3478627 ])
energies['gpu']['1440'] = np.array([  77987147,  77939565,  78052531 ])

print('\n')
energies_log = {}
for target in energies:
    energies_log[target] = {}
    energies_log[target]['avg'] = []
    for degree in energies[target]:
        values = energies[target][degree]
        avg    = np.mean(values)
        dev    = (np.max(values)-np.min(values))/avg
        print('For ' + str(target) + ' and degree=' + str(degree) + '...')
        print('  average energy usage is ' + str(round(avg,accuracy)) + ' s')
        print('  (max-min)/avg * 100 is ' + str(round(dev*100.0,accuracy)) + '.\n')
        energies_log[target]['avg'].append(np.log10(avg))


runtimes = {}
runtimes['cpu'] = {}
runtimes['gpu'] = {}

runtimes['cpu'][ 200] = np.array([    109,    109,    108 ])
runtimes['cpu'][ 300] = np.array([    340,    338,    337 ])
runtimes['cpu'][ 720] = np.array([   4397,   4366,   4377 ])
runtimes['cpu'][1440] = np.array([  40956,  40626,  40873 ])

runtimes['gpu'][ 200] = np.array([     23,     24,     23 ])
runtimes['gpu'][ 300] = np.array([     56,     56,     57 ])
runtimes['gpu'][ 720] = np.array([    966,    965,    964 ])
runtimes['gpu'][1440] = np.array([  11152,  11143,  11196 ])

print('\n')
runtimes_log = {}
for target in runtimes:
    runtimes_log[target] = {}
    runtimes_log[target]['avg'] = []
    for degree in runtimes[target]:
        values = runtimes[target][degree]
        avg    = np.mean(values)
        dev    = (np.max(values)-np.min(values))/avg
        print('For ' + str(target) + ' and degree=' + str(degree) + '...')
        print('  average runtime is ' + str(round(avg,accuracy)) + ' s')
        print('  (max-min)/avg * 100 is ' + str(round(dev*100.0,accuracy)) + '.\n')
        runtimes_log[target]['avg'].append(np.log10(avg)) 


shdeg  = [200, 300, 720, 1440]
shdeg_100 = [int(shd/100) if shd/100 % 1 == 0 else shd/100 for shd in shdeg]


fig = plt.gcf()


colour_cpu = '#66c2a5'
colour_gpu = '#fc8d62'

runtimes_log_diff = np.array(runtimes_log['cpu']['avg']) - np.array(runtimes_log['gpu']['avg'])
print('\nruntimes_log_diff: ' + str(runtimes_log_diff))

plt.plot(shdeg_100, runtimes_log['cpu']['avg'], linestyle='-', marker='o', color=colour_cpu, label='CPU')
plt.plot(shdeg_100, runtimes_log['gpu']['avg'], linestyle='-', marker='*', markersize=9, color=colour_gpu, label='GPU')

plt.xlabel('Degree / 100')
plt.xticks(shdeg_100, shdeg_100)
plt.ylabel('$Log_{10}(t)$')
plt.title(r'WMAM 5.0 Runtime on Cirrus')
plt.legend(loc="upper left")
plt.show()

fig.savefig('wmam-cpu-gpu-runtime.eps', format='eps', dpi=1000)
plt.show()
plt.clf()
fig = plt.gcf()


energies_log_diff = np.array(energies_log['cpu']['avg']) - np.array(energies_log['gpu']['avg'])
print('\nenergies_log_diff: ' + str(energies_log_diff))

energies_diff = np.array(energies_log['cpu']['avg']) - np.array(energies_log['gpu']['avg'])

plt.plot(shdeg_100, energies_log['cpu']['avg'], linestyle='-', marker='o', color=colour_cpu, label='CPU')
plt.plot(shdeg_100, energies_log['gpu']['avg'], linestyle='-', marker='*', markersize=9, color=colour_gpu, label='GPU')

plt.xlabel('Degree / 100')
plt.xticks(shdeg_100, shdeg_100)
plt.ylabel('$Log_{10}(E)$')
plt.title(r'WMAM 5.0 Energy Usage on Cirrus')
plt.legend(loc="upper left")
plt.show()

fig.savefig('wmam-cpu-gpu-energy.eps', format='eps', dpi=1000)
plt.show()
plt.clf()
fig = plt.gcf()
