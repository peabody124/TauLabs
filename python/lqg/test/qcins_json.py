#!/usr/bin/env python -B

# call with
#  python -m lqg.test.qcins_json PATH/TO/JSON

import argparse

parser = argparse.ArgumentParser(description='Process JSON log file and run through QC INS.')
parser.add_argument("file",
                    help  = "file to process")
# Parse the command-line.
args = parser.parse_args()
filename = args.file

import json
import numpy as np
import matplotlib.pyplot as plt
from lqg.test.qc_ins import init, advance, correct, configure
from ins.quaternions import quat_rpy

with open(filename) as infile:
    log = json.load(infile)

STEPS = len(log)

def get_field(fn):
    l = [(l['t'], l[fn]) for l in log if l[fn] is not None]
    t = [t for t,_ in l]
    rpy = np.array([rpy for _,rpy in l])
    return t,rpy

_,u = get_field('u')
beta_throttle = 1./np.median(u[:,-1],axis=0)

dT = 0.001
NUMX = 19

history = np.zeros((STEPS,NUMX))
rpy_history = np.zeros((STEPS,3))
history_t = np.zeros((STEPS,))

init()
configure(init_thrust=beta_throttle)
configure(init_bias=np.array([-0.1, 0.04, -0.2]))
#configure(sensor_noise=np.array([100., # altitude
#    1e6, 1e6, 1e6,                     # accel
#    1e3, 1e3, 1e3,                     # gyro
#    1e3, 1e3]),                        # mag
#    process_noise=np.array([
#        1e-3, 1e-3, 1e-3, 1e-4,  # position and velocity
#        1e-3, 1e-3, 1e-3, 1e-3,  # quaternion
#        1e-1, 1e-1, 1e-1,        # rotation rate
#        1e-3, 1e-3, 1e-3, 1e-1   # torques including throttle
#        ]),
#    mu=1)

for i in range(STEPS):

    params = dict()

    if False and log[i]['mag'] is not None:
        params['mag'] = np.array(log[i]['mag']) / np.sqrt(log[i]['mag'][0]**2 + log[i]['mag'][1] ** 2 + log[i]['mag'][2] ** 2)
        if i < 100:
            print(params)

    if False and log[i]['baro'] is not None:
        params['baro'] = -log[i]['baro']

    if i % 100 == 0: # fake baro data
        params['baro'] = 0.0

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

    if np.any(np.isnan(state)):
        print("NAN encountered")
        break

t_rpy,rpy = get_field('rpy')
t_gyros,gyros = get_field('gyros')
t_u,u = get_field('u')

#%matplotlib
#plt.figure(1)
ax1 = plt.subplot(3,4,1)
plt.plot(t_rpy, rpy[:,0], history_t, rpy_history[:,0])
plt.title('Roll Angle')
plt.subplot(3,4,2, sharex=ax1)
plt.plot(t_gyros, gyros[:,0], history_t, history[:,8] * 180.0 / np.pi)
plt.title('Roll Rate')
plt.subplot(3,4,3, sharex=ax1)
plt.plot(t_u, u[:,0], history_t, history[:,11])
plt.title('Roll Torque')
plt.subplot(3,4,4, sharex=ax1)
plt.plot(history_t, history[:,15])
plt.title('Roll Bias')

plt.subplot(3,4,5, sharex=ax1)
plt.plot(t_rpy, rpy[:,1], history_t, rpy_history[:,1])
plt.title('Pitch Angle')
plt.subplot(3,4,6, sharex=ax1)
plt.plot(t_gyros, gyros[:,1], history_t, history[:,9] * 180.0 / np.pi)
plt.title('Pitch Rate')
plt.subplot(3,4,7, sharex=ax1)
plt.plot(t_u, u[:,1], history_t, history[:,12])
plt.title('Pitch Torque')
plt.subplot(3,4,8, sharex=ax1)
plt.plot(history_t, history[:,16])
plt.title('Pitch Bias')

plt.subplot(3,4,9, sharex=ax1)
plt.plot(t_rpy, rpy[:,2], history_t, rpy_history[:,2])
plt.title('Yaw Angle')
plt.subplot(3,4,10, sharex=ax1)
plt.plot(t_gyros, gyros[:,2], history_t, history[:,10] * 180.0 / np.pi)
plt.subplot(3,4,11, sharex=ax1)
plt.title('Yaw Rate')
plt.plot(t_u, u[:,2], history_t, history[:,13])
plt.title('Yaw Torque')
plt.subplot(3,4,12, sharex=ax1)
plt.plot(history_t, history[:,17])
plt.title('Yaw Bias')

#plt.xlim(60,80)

t_vel, gps_vel = get_field('vel')
if len(gps_vel) > 0:
    plt.figure(2)
    #h_vel = np.sqrt(gps_vel[:,0]**2 + gps_vel[:,1]**2)
    #h_vel_ins = np.sqrt(history[:,1]**2 + history[:,2]**2)
    #plt.plot(t_vel, h_vel, history_t, h_vel_ins)
    ax1 = plt.subplot(2,1,1)
    plt.plot(t_vel, gps_vel[:,0], history_t, history[:,1])
    plt.title('North')
    ax1 = plt.subplot(2,1,2,sharex=ax1)
    plt.plot(t_vel, gps_vel[:,1], history_t, history[:,2])
    plt.title('East')


