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
import numpy.random as rand


def route_time(routes, route_c, t):
    """
    
    Parameters
    ----------
    routes : DataFrame
        Dataframe containing route information. Extracted from 
        Model_Parameters.xlsx.
        
    route_c : String (Numeric)
        The route number over which bus needs to be deployed.
    
    t : Float
        The clock time at deployment of the bus.

    Returns
    -------
    list(time_ar): List
        List object containing the times of arrival and service at 
        each stop of the route.
        
    list(stops): List
        List object containing the attribute of the travel/service times.
        1 denotes Travel time between stops and 0 denotes the service 
        time at a particular stop.
        
    """
    time_ar = []
    if route_c.isnumeric():
        search_r = 'route_' + route_c
    else:
        search_r = route_c
    sm_data = routes[routes.columns[routes.columns.str.match(search_r)]]
    sm_data.dropna(inplace=True)
    for i in range(sm_data.shape[0]):
        if sm_data[search_r+'_index'][i]==0:
            t_arr = rand.gamma(sm_data[search_r+'_mean'][i], 
                                 sm_data[search_r+'_std'][i])         # change distributions
            t = t + t_arr
        elif sm_data[search_r+'_index'][i]==1:
            t_ser = rand.gamma(sm_data[search_r+'_mean'][i], 
                                 sm_data[search_r+'_std'][i])         # change distributions
            t = t + t_ser
        time_ar.append(round(t, 2))
    stops = np.array(sm_data[search_r+'_index'])
    return list(time_ar), list(stops)


class bus:
    def __init__(self, charge, charge_std):          # Class initialization
        self.charge = round(rand.normal(charge, charge_std), 3)
        self.state = -1          # deployed = 1, refill = 0, standstill = -1
        self.route = None        # route in string, eg. '1', 'refill', 'recharge' 
        self.time_arr = list()   # array containing travel or stop service times for the bus
        self.event_arr = list()  # array to denote travel or stop service state
        
    def assign_route(self, routes, route_c, t):        # assigning a specific route to the bus
        self.time_arr, self.event_arr = route_time(routes, route_c, t)
        self.route = route_c
        if (route_c=='refill')or(route_c=='recharge'):      # deployed = 1, refill = 0, standstill = -1
            self.state = 0
        elif route_c.isnumeric():
            self.state = 1
    
    def next_t(self):            # passing the next event time
        if len(self.time_arr)==0:
            return np.inf
        else:
            return self.time_arr[0]
    
    def next_e(self):            # passing the next event type
        if len(self.event_arr)==0:
            return np.inf
        else:
            return self.event_arr[0]
        
    def last_t(self):            # passing the last event time for a route
        if len(self.event_arr)==0:
            return np.inf
        else:
            return self.time_arr[len(self.event_arr)-1]
    
    def info(self):
        bus_dict ={'charge': self.charge,
                   'state' : self.state,
                   'route' : self.route,
                   'event' : self.event_arr[0]}    # all details passed as a dictionary for table insertion
        return bus_dict


def gen_demands(routes, T):
    dem = pd.read_excel('Model_Parameters.xlsx', 'Demands')
    t_arr = np.array(range(0, T, 60))
    demand_r, demand_t, demand_c = [], [], []
    for row in dem.itertuples(index=False):
        if row[0]<=routes:
            # denerating demands
            route, a, b, c, d, charge = np.str_(row[0]), row[1], row[2], row[3], row[4], row[5]
            demand = np.ceil(a*np.sin((t_arr+c)/d) + b)             # a distribution can also be used
            # generating times wrt demands
            for i, t in enumerate(t_arr):
                if demand[i]>0:
                    for t_d in range(0, 60, int(60/demand[i])):
                        demand_t.append(round(t + t_d, 2))
                        demand_r.append(route)
                        demand_c.append(charge)

    demand_r.append(None)
    demand_c.append(np.inf)
    demand_t.append(np.inf)
    all_dt = pd.DataFrame(columns=['routes', 'times', 'charge'])
    all_dt['routes'] = demand_r
    all_dt['times'] = demand_t
    all_dt['charge'] = demand_c
    all_dt = all_dt.sort_values(['times','routes'], ignore_index=True)
    demand_r = np.array(all_dt['routes'])
    demand_c = np.array(all_dt['charge'])
    demand_t = np.array(all_dt['times'])
    return list(demand_t), list(demand_r), list(demand_c)


# +
def next_bus_e(buses):
    min_t = np.inf
    index = None
    for i in range(len(buses)):
        if (buses[i].next_t()<min_t):
            min_t = buses[i].next_t()
            index = i
    return min_t, buses[i].next_e(), index

def available_bus(buses, dem_charge):
    b_charges = [buses[i].charge for i in range(len(buses)) if buses[i].state==-1]
    b_index = [i for i in range(len(buses)) if buses[i].state==-1]
    index = -1
    if len(b_charges)>0:
        if (max(b_charges)>dem_charge):
            for i in range(len(b_charges)):
                if max(b_charges)==b_charges[i]:
                    index = b_index[i]
    return index

def unavailable_bus(buses, min_charge):
    b_charges = [buses[i].charge for i in range(len(buses)) if buses[i].state==-1]
    b_index = [i for i in range(len(buses)) if buses[i].state==-1]
    index = -1
    if (len(b_charges)>0)and(min_charge!=np.inf):
        if (min(b_charges)<min_charge):
            for i in range(len(b_charges)):
                if min(b_charges)==b_charges[i]:
                    index = b_index[i]
    return index

def buses_status(buses):
    n_dep = sum([1 for i in buses if i.state==1])
    n_ref = sum([1 for i in buses if i.state==0])
    n_stds = sum([1 for i in buses if i.state==-1])
    return n_dep, n_ref, n_stds


# +
def SS_update(SS_table, t, buses, dct, bus_e, t_updt, index=-1, dem_ct=0, dem_at=0, dem_c=0):
    ss_d = {}
    ss_d['Time'] = t_updt
    ss_arr = [len(buses), dct]
    ss_dep, ss_re, ss_ss = buses_status(buses)
    ss_arr.append(ss_dep)
    ss_arr.append(ss_re)
    ss_arr.append(ss_ss)
    ss_d['System_State'] = ss_arr
    if index!=-1:
        ss_d['Bus'] = index + 1
        ss_d['Charge'] = buses[index].charge
        ss_d['Route'] = buses[index].route
        ss_d['State'] = buses[index].state
    else:
        ss_d['Bus'] = np.nan
        ss_d['Charge'] = np.nan
        ss_d['Route'] = None
        ss_d['State'] = np.nan
    ss_d['Event'] = bus_e
    ss_d['Process_Time'] = t_updt - t    
    if dem_c==0:
        ss_d['Demand_Current'] = np.nan
        ss_d['Demand_Actual'] = np.nan
        ss_d['Demand_Charge'] = np.nan
    else:
        ss_d['Demand_Current'] = dem_ct
        ss_d['Demand_Actual'] = dem_at
        ss_d['Demand_Charge'] = dem_c
    SS_table = SS_table.append(ss_d, ignore_index=True)
    return SS_table

def BD_update(BD_table, t_updt, dem_ct, dem_at, dem_c, buses, index):
    bd_d = {}
    bd_d['Time'] = t_updt
    bd_d['Demand_Current'] = dem_ct
    bd_d['Demand_Actual'] = dem_at
    bd_d['Demand_Charge'] = dem_c
    bd_d['Bus'] = index + 1
    bd_d['Charge'] = buses[index].charge
    bd_d['Route'] = buses[index].route
    bd_d['State'] = buses[index].state
    bd_d['Event_Array'] = np.array(buses[index].event_arr)
    bd_d['Time_Array'] = np.array(buses[index].time_arr)
    BD_table = BD_table.append(bd_d, ignore_index=True)
    return BD_table
# -


