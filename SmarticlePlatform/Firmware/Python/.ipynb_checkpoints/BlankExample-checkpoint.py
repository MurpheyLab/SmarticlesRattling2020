# SmarticleSyncExample.py
import sys
sys.path.append('modules')

from SmarticleSwarm import *
from random import randint
import time
import numpy as np
import math

# TO DO: change to your port
PORT_NAME = '/dev/tty.usbserial-DN050I6Q'


def rx_callback(xbee_message):
    '''Simple call_back function to print received packets'''

    # print("From {}  {}".format(xbee_message.remote_device.get_node_id(),
    #                          xbee_message.data.decode("utf-8",errors="replace")))
    print(xbee_message.data.decode("utf-8",errors="backslashreplace"), end ='')


swarm = SmarticleSwarm(port=PORT_NAME, debug=1)
swarm.xb.add_rx_callback(rx_callback)
swarm.build_network(8)

swarm.set_mode(2)
L1 = [0,180,180,0]
R1 = [0,0,180,180]

L2 = [45,135,135,45]
R2 = [45,45,135,135]
swarm.gait_init([L1,R1], 450,gait_num = 0) #gaits, delay between poitns in ms; I wouldnt go faster than 200m
# swarm.gait_init([L2,R2],450, gait_num = 1) #gaits, delay between poitns in ms; I wouldnt go faster than 200ms
# swarm.init_sync_thread()
