# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.12.0
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

import pandas as pd
import numpy as  np
import string as str
import numpy.random as random


def route_time(routes, route_c, t):
    time_ar = []
    if route_c.isnumeric():
        search_r = 'route_' + route_c
    else:
        search_r = route_c
    sm_data = routes[routes.columns[routes.columns.str.match(search_r)]]
    sm_data.dropna(inplace=True)
    for i in range(sm_data.shape[0]):
        if sm_data[search_r+'_index'][i]==0:
            t_arr = random.gamma(sm_data[search_r+'_mean'][i], sm_data[search_r+'_std'][i])
            t = t + t_arr
        elif sm_data[search_r+'_index'][i]==1:
            t_ser = random.gamma(sm_data[search_r+'_mean'][i], sm_data[search_r+'_std'][i])
            t = t + t_ser
        time_ar.append(round(t, 4))
    stops = np.array(sm_data[search_r+'_index'])
    return time_ar, stops


class bus:
    def __init__(self):
        self.charge = random
        self.state = -1      # deployed = 1, refill = 0, standstill = -1
        self.route = None    # route in string, eg. '1', 'refill', 'recharge' 
        self.time_arr = []   # array containing travel or stop service times for the bus
        self.event_arr = []   # array to denote travel or stop service state
        
    def assign_route(self, route_t, route_c, t):
        self.time_arr, self.event_arr = route_time(route_t, route_c, t)
    
    def next_t(self):
        return self.time_arr[0]
    def next_e(self):
        return self.event_arr[0]
