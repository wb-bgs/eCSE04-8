#!/usr/bin/env python

'''
This Python script converts nvidia-smi power readings to energy consumed per GPU.

It assumes that a command of the following form has been run.

nvidia-smi --query-gpu=index,timestamp,power.draw --format=csv --loop=1 &> nvsmi-power-${SLURM_JOB_ID}.out &

The resulting 'nvsmi-power-${SLURM_JOB_ID}.out' file is then assumed to contain text in the following format.

index, timestamp, power.draw [W]
0, 2025/04/09 08:59:06.356, 63.86 W
1, 2025/04/09 08:59:06.360, 56.78 W
2, 2025/04/09 08:59:06.362, 56.33 W
3, 2025/04/09 08:59:06.363, 59.83 W
...
'''



import sys
import os
import datetime
import argparse
import glob

import numpy as np
import scipy


parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('-v',  '--version', action='version', version='%(prog)s 1.0')
parser.add_argument('-p',  '--power-readings-mask', type=str, help='The path mask for files containing the power readings read via an nvidia-smi query command.')
parser.add_argument('-a',  '--accuracy', type=int, default=3, help='The number of decimal places used for rounding.')

args = parser.parse_args()


NVSMI_GPUNUM_INDEX    = 0
NVSMI_TIMESTAMP_INDEX = 1
NVSMI_POWER_INDEX     = 2
NVSMI_FIELD_COUNT     = 3


totalEnergySum = 0.0

for power_readings_fn in glob.glob(args.power_readings_mask):

    timeStamps = {}
    power = {}

    print('\nReading ' + str(power_readings_fn) + '...\n')
    nvsmiData  = os.popen('cat ' + power_readings_fn).read()
    nvsmiLines = nvsmiData.split('\n')

    for ln in nvsmiLines:
        if ln[:5] == 'index':
            # line is header
            continue

        lnParts = ln.split(',')
        if len(lnParts) < NVSMI_FIELD_COUNT:
            continue
        
        gpuNum = int(lnParts[NVSMI_GPUNUM_INDEX])
        if gpuNum not in timeStamps:
            timeStamps[gpuNum] = []
            power[gpuNum] = []

        timeStamps[gpuNum].append(datetime.datetime.strptime(lnParts[NVSMI_TIMESTAMP_INDEX][1:], '%Y/%m/%d %H:%M:%S.%f'))
        power[gpuNum].append(float(lnParts[NVSMI_POWER_INDEX].split()[0]))


    gpus = len(timeStamps)

    times = {}
    for i in range(gpus):
        times[i] = [0.0]
        tmstmp0 = timeStamps[i][0]
        for tmstmp in timeStamps[i][1:]:
            tmdelta = tmstmp - tmstmp0
            times[i].append(tmdelta.total_seconds())

    energy = {}
    totalEnergy = 0.0
    for i in range(gpus):
        x = np.array(times[i])
        y = np.array(power[i])
        energy[i] = scipy.integrate.cumulative_trapezoid(y, x, initial=0)

        gpuEnergy = energy[i][-1]
        print('GPU ' + str(i) + ' consumed ' + str(round(gpuEnergy,args.accuracy)) + ' J.')
        totalEnergy += gpuEnergy

    print('\nTotal energy consumed by GPU node is ' + str(round(totalEnergy,args.accuracy)) +  ' J.')
    totalEnergySum += totalEnergy


print('\n\nTotal energy consumed by all GPU nodes is ' + str(round(totalEnergySum,args.accuracy)) +  ' J.')