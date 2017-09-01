#!/usr/bin/env python -B

"""
Extract from an Overo log the core updates that would be useful for 
fitting that data with new filters and store that as JSON for use
later.

Calling format:
 ./lqg/ins_fit_log.py -t -g <GITHASH> <FILENAME>
"""

from taulabs import telemetry, uavo
uavo_list = telemetry.get_telemetry_by_args()

accels = uavo_list.as_numpy_array(uavo.UAVO_Accels)
gyros = uavo_list.as_numpy_array(uavo.UAVO_Gyros)
mags = uavo_list.as_numpy_array(uavo.UAVO_Magnetometer)
baros = uavo_list.as_numpy_array(uavo.UAVO_BaroAltitude)
actuator = uavo_list.as_numpy_array(uavo.UAVO_ActuatorDesired)
attitude = uavo_list.as_numpy_array(uavo.UAVO_AttitudeActual)
gps_vel = uavo_list.as_numpy_array(uavo.UAVO_GPSVelocity)
gps_pos = uavo_list.as_numpy_array(uavo.UAVO_GPSVelocity)


from matplotlib.mlab import find

# find the time of the last gyro update before the last actuator update
# to ensure we always have valid data for running
t_start = actuator['time'][find(actuator['Throttle'] > 0)[0]]+.1
t_end = actuator['time'][find(actuator['Throttle'] > 0)[-1]]

g_idx_start = find(gyros['time'] > t_start)[0]
g_idx_end = find(gyros['time'] < t_end)[-1]
g_idx = np.arange(g_idx_start, g_idx_end)
STEPS = len(g_idx)
print("Length: {}".format(STEPS))

last_b_idx = -1
last_m_idx = -1
last_gps_vel_idx = -1
last_gps_pos_idx = -1

log = []

for i in range(STEPS):

    idx = g_idx[i] # look up the gyro index for this time step
    t = gyros['time'][idx]
    c_idx = find(actuator['time'] < t)[-1]  # find index for the most recent control before this gyro
    a_idx = find(accels['time'] < t)[-1]    # find index for the most recent accel before this gyro
    b_idx = find(baros['time'] < t)[-1]     # find index for the most recent baro before this gyro
    m_idx = find(mags['time'] < t)[-1]      # find index for the most recent mag before this gyro
    att_idx = find(attitude['time'] < t)[-1]
    
    g = [gyros['x'][i,0], gyros['y'][i,0], gyros['z'][i,0]],
    a = [accels['x'][a_idx,0], accels['y'][a_idx,0], accels['z'][a_idx,0]]
    u = [actuator['Roll'][c_idx,0], actuator['Pitch'][c_idx,0], actuator['Yaw'][c_idx,0], actuator['Throttle'][c_idx,0]]
    rpy = [attitude['Roll'][att_idx,0], attitude['Pitch'][att_idx,0], attitude['Yaw'][att_idx,0]]

    baro = None
    if b_idx > last_b_idx:
        baro = baros['Altitude'][b_idx,0]
        last_b_idx = b_idx

    mag = None
    if m_idx > last_m_idx:
        mag = [mags['x'][m_idx,0], mags['y'][m_idx,0], mags['z'][m_idx,0]]
        last_m_idx = m_idx

    gps_pos_val = None
    if gps_pos is not None and any(gps_pos['time'] < t):
        gps_pos_idx = find(gps_pos['time'] < t)[-1]
        if gps_pos_idx > last_gps_pos_idx:
            gps_pos_val = [gps_pos['North'][gps_pos_idx,0], gps_pos['East'][gps_pos_idx,0], gps_pos['Down'][gps_pos_idx,0]]
            last_gps_pos_idx = gps_pos_idx

    gps_vel_val = None
    if gps_vel is not None and any(gps_vel['time'] < t):
        gps_vel_idx = find(gps_vel['time'] < t)[-1]
        if gps_vel_idx > last_gps_vel_idx:
            gps_vel_val = [gps_vel['North'][gps_vel_idx,0], gps_vel['East'][gps_vel_idx,0], gps_vel['Down'][gps_vel_idx,0]]
            last_gps_vel_idx = gps_vel_idx


    log.append({'t': t, 'gyros': g, 'accels': a, 'baro': baro, 'mag': mag, 'pos': gps_pos_val, 'vel': gps_vel_val, 'rpy': rpy})

filename = "temp.JSON"

import io
try:
    to_unicode = unicode
except NameError:
    to_unicode = str

with io.open(filename, 'w', encoding='utf8') as outfile:
    str_ = json.dumps(log, separators=(',', ': '), ensure_ascii=False)
    outfile.write(to_unicode(str_))

#with open(filename) as infile:
#    log_check = json.load(infile)