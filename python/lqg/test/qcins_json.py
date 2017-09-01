#!/usr/bin/env python -B

import numpy as np
import matplotlib.pyplot as plt

import json
filename = "temp.JSON"
with open(filename) as infile:
    log = json.load(infile)

STEPS = len(log)

beta_throttle = 1.4 # 1.0/np.median(actuator['Throttle'])

dT = 0.001
NUMX = 15

history = np.zeros((STEPS,NUMX))
rpy_history = np.zeros((STEPS,3))
history_t = np.zeros((STEPS,))

from lqg.test.qc_ins import init, advance, correct, configure
from ins.quaternions import quat_rpy

init()
for i in range(STEPS):

    params = dict()

    if log[i]['mag'] is not None:
        params['mag'] = np.array(log[i]['mag'])

    if log[i]['baro'] is not None:
        params['baro'] = log[i]['baro']

    # TODO: pull out gps position and velocity
    
    g = np.array(log[i]['gyros']) * np.pi / 180.0
    a = np.array(log[i]['accels'])
    u = np.array(log[i]['u'])

    state = advance(u, dT)
    state = correct(g, a, **params)

    history[i,:] = state
    history_t[i] = log[i]['t']

    rpy_history[i,:] = quat_rpy(state[4:8])

    if i % 1000 == 0:
        print("{} / {}".format(i, STEPS))

def get_field(fn):
    l = [(l['t'], l[fn]) for l in log if l[fn] is not None]
    t = [t for t,_ in l]
    rpy = np.array([rpy for _,rpy in l])
    return t,rpy

t_rpy,rpy = get_field('rpy')
t_gyros,gyros = get_field('gyros')

ax1 = plt.subplot(3,2,1)
plt.plot(t_rpy, rpy[:,0], history_t, rpy_history[:,0])
plt.title('Roll')
plt.subplot(3,2,2, sharex=ax1)
plt.plot(t_gyros, gyros[:,0], history_t, history[:,8] * 180.0 / np.pi)
plt.subplot(3,2,3, sharex=ax1)
plt.plot(t_rpy, rpy[:,1], history_t, rpy_history[:,1])
plt.title('Pitch')
plt.subplot(3,2,4, sharex=ax1)
plt.plot(t_gyros, gyros[:,1], history_t, history[:,9] * 180.0 / np.pi)
plt.subplot(3,2,5, sharex=ax1)
plt.plot(t_rpy, rpy[:,2], history_t, rpy_history[:,2])
plt.title('Yaw')
plt.subplot(3,2,6, sharex=ax1)
plt.plot(t_gyros, gyros[:,2], history_t, history[:,10] * 180.0 / np.pi)

#plt.figure()
#gps_vel_ne = sqrt(square(gps_vel['North']) + square(gps_vel['East']))
#vel = np.sqrt(history[:,1]**2 + history[:,2]**2)
#plot(gps_vel['time'], gps_vel_ne, history_t, vel)

