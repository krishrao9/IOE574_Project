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
import matplotlib.pyplot as plt
import scipy.stats as st
import statsmodels.api as sm


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
    def __init__(self, charge, charge_std):            # Class initialization
        self.charge = round(rand.normal(charge, charge_std), 3)
        self.state = -1                                # deployed = 1, refill = 0, standstill = -1
        self.route = None                              # route in string, eg. '1', 'refill', 'recharge' 
        self.time_arr = list()                         # array containing travel or stop service times for the bus
        self.event_arr = list()                        # array to denote travel or stop service state
        
    def assign_route(self, routes, route_c, t):        # assigning a specific route to the bus
        self.time_arr, self.event_arr = route_time(routes, route_c, t)
        self.route = route_c
        if (route_c=='refill')or(route_c=='recharge'): # deployed = 1, refill = 0, standstill = -1
            self.state = 0
        elif route_c.isnumeric():
            self.state = 1
    
    def assign_route_varred(self, e_array, t_array, route_c, t):  # assigning a specific route to the bus
        self.time_arr, self.event_arr = list(t_array), list(e_array)
        self.route = route_c
        if (route_c=='refill')or(route_c=='recharge'):  # deployed = 1, refill = 0, standstill = -1
            self.state = 0
        elif route_c.isnumeric():
            self.state = 1
    
    def next_t(self):                       # passing the next event time
        if len(self.time_arr)==0:
            return np.inf
        else:
            return self.time_arr[0]
    
    def next_e(self):                       # passing the next event type
        if len(self.event_arr)==0:
            return np.inf
        else:
            return self.event_arr[0]
        
    def last_t(self):                       # passing the last event time for a route
        if len(self.event_arr)==0:
            return np.inf
        else:
            return self.time_arr[len(self.event_arr)-1]
    
    def info(self):                         # all details passed as a dictionary for table insertion
        bus_dict = {'charge': self.charge,
                    'state' : self.state,
                    'route' : self.route,
                    'event' : self.event_arr[0]}
        return bus_dict


def gen_demands(routes, T):
    dem = pd.read_excel('Model_Parameters.xlsx', 'Demands')
    t_arr = np.array(range(0, T, 60))
    demand_r, demand_t, demand_c = [], [], []
    for row in dem.itertuples(index=False):
        if row[0]<=routes:
        # generating demands
            route, a, b, c, d, charge = np.str_(row[0]), row[1], row[2], row[3], row[4], row[5]
            demand = np.ceil(a*np.sin((t_arr+c)/d) + b)        # a distribution can also be used
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
def fleet_simulation(t, T, routes, buses, refuel, running_consumption, service_consumption, 
                     refuel_stations, refuel_consumption, conversion_factor,
                     demand_at, demand_ct, demand_r, demand_c, dct_flag, SS_cols, BD_cols):
    ss_table = pd.DataFrame(columns=SS_cols)
    bd_table = pd.DataFrame(columns=BD_cols)

    # initial SS update
    time_check = np.inf
    ss_table = SS_update(ss_table, t, buses, dct_flag, np.nan, np.nan)
    
    # previous times for fleet
    prev_time = [np.nan]*len(buses)
    
    
    t_check = []
    # Simulate! Simulate! Simulate!        
    while (t<T)or(time_check!=np.inf):
        #print('---\nNew Event')
        time_check = min(next_bus_e(buses)[0], demand_ct[0])
        
        # -----
        # Priority One: Updating demands
        if ((demand_ct[0]!=np.inf) and
            (demand_ct[0]==time_check)):
            new_demand = [1 for i in range(len(demand_ct)) if demand_ct[i]==time_check]
            dct_flag = sum(new_demand)
            # testing
            t = time_check
        
        #----------
        # Main Switch Statements
        bus_chk = available_bus(buses, demand_c[0])
        
        #-----
        # Case 1: Sending low-fuel buses to refuel
        refuel_index = unavailable_bus(buses, min(demand_c))
        if ((time_check!=np.inf) and ((t<T) or (dct_flag>0)) and
            (refuel_index!=-1) and
            (buses_status(buses)[1]<refuel_stations)): 
            # and(refuel_index!=np.inf)):- might create complications
            dem_ct, dem_at, dem_c = np.nan, np.nan, np.nan
            #print('\tSending Bus for Refuel')
            msg = 'Sending Bus for Refuel'
            buses[refuel_index].assign_route(routes, refuel, t)
            prev_time[refuel_index] = t
            bus_e = 0
            t_updt = t

            ss_table = SS_update(ss_table, t, buses, dct_flag, bus_e, t_updt, refuel_index)
            bd_table = BD_update(bd_table, t_updt, dem_ct, dem_at, dem_c, buses, refuel_index)

            t = t_updt
            refuel_index = -1


        #-----
        # Case 2: Checking bus availability and deploying buses
        elif ((t<T) and
              (demand_ct[0]==time_check) and
              (demand_ct[0]!=np.inf) and
              (bus_chk!=-1)):

            #print('\tDeploying Bus')
            msg = 'Deploying Bus'
            dem_ct, dem_at = demand_ct.pop(0), demand_at.pop(0)
            dem_c, dem_r = demand_c.pop(0), demand_r.pop(0)
            index = available_bus(buses, dem_c)
            t_updt = dem_ct
            buses[index].assign_route(routes, dem_r, t_updt)
            prev_time[index] = t_updt
            bus_e = 0

            ss_table = SS_update(ss_table, t, buses, dct_flag, bus_e, t_updt, index, dem_ct, dem_at, dem_c)
            bd_table = BD_update(bd_table, t_updt, dem_ct, dem_at, dem_c, buses, index)
            dct_flag -= 1
            t = t_updt

            dem_ct, dem_at, dem_c, dem_r = np.nan, np.nan, np.nan, np.nan
            t_updt = np.nan
            index = np.nan


        #-----
        # Case 3: Checking next bus event and updating SS
        elif (((t<T) and
              (next_bus_e(buses)[0]==time_check) and 
              (next_bus_e(buses)[0]!=np.inf) and 
              ((buses_status(buses)[0]>0) or (buses_status(buses)[1]>0))) or
              (((t<T) and
              (demand_ct[0]==time_check) and 
              (demand_ct[0]!=np.inf) and 
              (bus_chk==-1)))):

            #print('\tUpdating next bus event')
            msg = 'Updating next bus event'
            index = next_bus_e(buses)[2]
            t_updt = buses[index].time_arr.pop(0)
            bus_e = buses[index].event_arr.pop(0)
            diff_t = prev_time[index]
            
            if (bus_e==1):
                mul_c = running_consumption
            elif (bus_e==0)and(buses[index].state==1):
                mul_c = service_consumption
            elif (bus_e==0)and(buses[index].state==0):
                mul_c = refuel_consumption
            buses[index].charge = round(buses[index].charge - (t_updt - diff_t)*mul_c*conversion_factor, 3)

            ss_table = SS_update(ss_table, diff_t, buses, dct_flag, bus_e, t_updt, index)
            t = t_updt
            prev_time[index] = t
            
            # if its the last event for the bus, update bus parameters to standstill
            if buses[index].next_e()==np.inf:
                buses[index].state = -1
                buses[index].route = None
                buses[index].time_arr = list()
                buses[index].event_arr = list()
                prev_time[index] = np.nan
            index = np.nan
            t_updt = np.nan
            bus_e = np.nan
        
              
        #-----
        # Case 4: All demands are completed within the timeframe 
        elif ((t<T) and
              (time_check==np.inf)):
            #print('\tJumping to EOD')
            msg = 'Jumping to EOD'
            t_updt = T
            ss_table = SS_update(ss_table, t, buses, dct_flag, bus_e, t_updt)
            t = t_updt    


        #-----
        # Case 5: Demands still exist beyond timeframe and need to be met
        elif ((t>=T) and
              (demand_ct[0]==time_check) and
              (demand_ct[0]!=np.inf) and
              (bus_chk!=-1)):
            #print('\tDeploying Bus beyond timeframe')
            msg = 'Deploying Bus beyond timeframe'
            dem_ct, dem_at = demand_ct.pop(0), demand_at.pop(0)
            dem_c, dem_r = demand_c.pop(0), demand_r.pop(0)
            index = available_bus(buses, dem_c)
            t_updt = dem_ct
            buses[index].assign_route(routes, dem_r, t_updt)
            prev_time[index] = t_updt
            bus_e = 0
            ss_table = SS_update(ss_table, t, buses, dct_flag, bus_e, t_updt, index, dem_ct, dem_at, dem_c)
            bd_table = BD_update(bd_table, t_updt, dem_ct, dem_at, dem_c, buses, index)
            dct_flag -= 1
            t = t_updt
            dem_ct, dem_at, dem_c, dem_r = np.nan, np.nan, np.nan, np.nan
            t_updt = np.nan
            index = np.nan


        #-----
        # Case 6: Going beyond timeframe, checking deployed buses and updating SS
        elif (((t>=T) and
              (next_bus_e(buses)[0]==time_check) and
              (next_bus_e(buses)[0]!=np.inf)) or
              ((t>=T) and
              (demand_ct[0]==time_check) and 
              (demand_ct[0]!=np.inf) and 
              (bus_chk==-1))):
            #print('\tUpdating next bus event beyond timeframe')
            msg = 'Updating next bus event beyond timeframe'
            index = next_bus_e(buses)[2]
            t_updt = buses[index].time_arr.pop(0)
            bus_e = buses[index].event_arr.pop(0)
            diff_t = prev_time[index]
            if (bus_e==1):
                mul_c = running_consumption
            elif (bus_e==0)and(buses[index].state==1):
                mul_c = service_consumption
            elif (bus_e==0)and(buses[index].state==0):                       # different for 'refill' and 'recharge' 
                mul_c = refuel_consumption                                   # difference from 90 for refill
            buses[index].charge = buses[index].charge - (t_updt - diff_t)*mul_c*conversion_factor
            ss_table = SS_update(ss_table, diff_t, buses, dct_flag, bus_e, t_updt, index)
            t = t_updt
            prev_time[index] = t
            
            # if its the last event for the bus, update bus parameters to standstill
            if buses[index].next_e()==np.inf:
                buses[index].state = -1
                buses[index].route = None
                buses[index].time_arr = list()
                buses[index].event_arr = list()
                prev_time[index] = np.nan
            index = np.nan
            t_updt = np.nan
            bus_e = np.nan

        
        #-----
        # Case 7: All demands are completed outside the timeframe 
        elif ((t>=T) and
              (time_check==np.inf)):
            #print('\tJumping to EOD')
            msg = 'Jumping to EOD'
            t_updt = np.inf        
            ss_table = SS_update(ss_table, t, buses, dct_flag, bus_e, t_updt)
            t = t_updt    

        # updating the current demand - if demand exists it'll automatically be updated
        demand_ct = list(np.sort(demand_ct))
        for k in range(dct_flag):
            demand_ct[k] = t
        
        # code to check if the program is running in loops
        if len(t_check)<5:
            t_check.append(t)
        else:
            t_check.pop(0)
            t_check.append(t)
        if (len(t_check)==5)and(t_check[0]==max(t_check)):
            print('Stuck in a loop!')
            print(msg)
            ss_table.to_parquet('ss_check.parquet')
            bd_table.to_parquet('bd_check.parquet')
            print('\tTime -', t)
            print('\tTime Check -', time_check)
            print('\tDemands -', demand_ct[:6], demand_at[:6], demand_r[:6], demand_c[:6])
            print('\tdct_flag -', dct_flag)
        
    return ss_table, bd_table


def fleet_simulation_varred(t, T, routes, event_array, time_array, init_time, buses, refuel, running_consumption, 
                     service_consumption, refuel_stations, refuel_consumption, conversion_factor,
                     demand_at, demand_ct, demand_r, demand_c, dct_flag, SS_cols, BD_cols):
    
    ss_table = pd.DataFrame(columns=SS_cols)
    bd_table = pd.DataFrame(columns=BD_cols)

    # Initial SS update
    time_check = np.inf
    ss_table = SS_update(ss_table, t, buses, dct_flag, np.nan, np.nan)
    
    # previous times for fleet
    prev_time = [np.nan]*len(buses)
    
    # Simulate! Simulate! Simulate!        
    while (t<T)or(time_check!=np.inf):
        #print('---\nNew Event')
        time_check = min(next_bus_e(buses)[0], demand_ct[0])
        
        # -----
        # Priority One: Updating demands
        if ((demand_ct[0]!=np.inf) and
            (demand_ct[0]==time_check)):
            new_demand = [1 for i in range(len(demand_ct)) if demand_ct[i]==time_check]
            dct_flag = sum(new_demand)
            t = time_check
        
        #----------
        # Main Simulation
        bus_chk = available_bus(buses, demand_c[0])

        #-----
        # Case 1: Sending low-fuel buses to refuel
        refuel_index = unavailable_bus(buses, min(demand_c))
        if (((t<T) or(dct_flag>0)) and  # or (time_check!=np.inf)
            (refuel_index!=-1) and
            (buses_status(buses)[1]<refuel_stations)): 
            dem_ct, dem_at, dem_c = np.nan, np.nan, np.nan
            #print('\tSending Bus for Refuel')
            buses[refuel_index].assign_route(routes, refuel, t)
            prev_time[refuel_index] = t
            bus_e = 0
            t_updt = t

            ss_table = SS_update(ss_table, t, buses, dct_flag, bus_e, t_updt, refuel_index)
            bd_table = BD_update(bd_table, t_updt, dem_ct, dem_at, dem_c, buses, refuel_index)

            t = t_updt
            refuel_index = -1


        #-----
        # Case 2: Checking bus availability and deploying buses
        elif ((t<T) and
              (demand_ct[0]==time_check) and
              (demand_ct[0]!=np.inf) and
              (bus_chk!=-1)):

            #print('\tDeploying Bus')
            dem_ct, dem_at = demand_ct.pop(0), demand_at.pop(0)
            dem_c, dem_r = demand_c.pop(0), demand_r.pop(0)
            index = available_bus(buses, dem_c)
            t_updt = dem_ct
            be_arr = event_array.pop(0)
            bt_arr = np.array(time_array.pop(0))
            bit_arr = dem_ct - np.array(init_time.pop(0))
            buses[index].assign_route_varred(be_arr, bt_arr+bit_arr, dem_r, t_updt)
            prev_time[index] = t_updt
            bus_e = 0

            ss_table = SS_update(ss_table, t, buses, dct_flag, bus_e, t_updt, index, dem_ct, dem_at, dem_c)
            bd_table = BD_update(bd_table, t_updt, dem_ct, dem_at, dem_c, buses, index)
            dct_flag -= 1
            t = t_updt

            dem_ct, dem_at, dem_c, dem_r = np.nan, np.nan, np.nan, np.nan
            t_updt = np.nan
            index = np.nan


        #-----
        # Case 3: Checking next bus event and updating SS
        elif (((t<T) and
              (next_bus_e(buses)[0]==time_check) and 
              (next_bus_e(buses)[0]!=np.inf) and 
              ((buses_status(buses)[0]>0) or (buses_status(buses)[1]>0))) or
              (((t<T) and
              (demand_ct[0]==time_check) and 
              (demand_ct[0]!=np.inf) and 
              (bus_chk==-1)))):

            #print('\tUpdating next bus event')
            index = next_bus_e(buses)[2]
            t_updt = buses[index].time_arr.pop(0)
            bus_e = buses[index].event_arr.pop(0)
            diff_t = prev_time[index]
            
            if (bus_e==1):
                mul_c = running_consumption
            elif (bus_e==0)and(buses[index].state==1):
                mul_c = service_consumption
            elif (bus_e==0)and(buses[index].state==0):                       # different for 'refill' and 'recharge' 
                mul_c = refuel_consumption                                   # difference from 90 for refill
            buses[index].charge = round(buses[index].charge - (t_updt - diff_t)*mul_c*conversion_factor, 3)

            ss_table = SS_update(ss_table, diff_t, buses, dct_flag, bus_e, t_updt, index)
            t = t_updt
            prev_time[index] = t
            
            # if its the last event for the bus, update bus parameters to standstill
            if buses[index].next_e()==np.inf:
                buses[index].state = -1
                buses[index].route = None
                buses[index].time_arr = list()
                buses[index].event_arr = list()
                prev_time[index] = np.nan
            index = np.nan
            t_updt = np.nan
            bus_e = np.nan
        
        
        #-----
        # Case 4: All demands are completed within the timeframe 
        elif ((t<T) and
              (time_check==np.inf)):

            #print('\tJumping to EOD')
            t_updt = T

            ss_table = SS_update(ss_table, t, buses, dct_flag, bus_e, t_updt)
            t = t_updt    


        #-----
        # Case 5: Demands still exist beyond timeframe and need to be met
        elif ((t>=T) and
              (demand_ct[0]==time_check) and
              (demand_ct[0]!=np.inf) and
              (bus_chk!=-1)):

            #print('\tDeploying Bus beyond timeframe')
            dem_ct, dem_at = demand_ct.pop(0), demand_at.pop(0)
            dem_c, dem_r = demand_c.pop(0), demand_r.pop(0)
            index = available_bus(buses, dem_c)
            t_updt = dem_ct
            be_arr = event_array.pop(0)
            bt_arr = np.array(time_array.pop(0))
            bit_arr = dem_ct - np.array(init_time.pop(0))
            buses[index].assign_route_varred(be_arr, bt_arr+bit_arr, dem_r, t_updt)
            prev_time[index] = t_updt
            bus_e = 0

            ss_table = SS_update(ss_table, t, buses, dct_flag, bus_e, t_updt, index, dem_ct, dem_at, dem_c)
            bd_table = BD_update(bd_table, t_updt, dem_ct, dem_at, dem_c, buses, index)
            dct_flag -= 1
            t = t_updt

            dem_ct, dem_at, dem_c, dem_r = np.nan, np.nan, np.nan, np.nan
            t_updt = np.nan
            index = np.nan


        #-----
        # Case 6: Going beyond timeframe, checking deployed buses and updating SS
        elif (((t>=T) and
              (next_bus_e(buses)[0]==time_check) and
              (next_bus_e(buses)[0]!=np.inf)) or
              ((t>=T) and
              (demand_ct[0]==time_check) and 
              (demand_ct[0]!=np.inf) and 
              (bus_chk==-1))):
            
            #print('\tUpdating next bus event beyond timeframe')
            index = next_bus_e(buses)[2]
            t_updt = buses[index].time_arr.pop(0)
            bus_e = buses[index].event_arr.pop(0)
            diff_t = prev_time[index]
            
            if (bus_e==1):
                mul_c = running_consumption
            elif (bus_e==0)and(buses[index].state==1):
                mul_c = service_consumption
            elif (bus_e==0)and(buses[index].state==0):                       # different for 'refill' and 'recharge' 
                mul_c = refuel_consumption                                   # difference from 90 for refill
            buses[index].charge = buses[index].charge - (t_updt - diff_t)*mul_c*conversion_factor

            ss_table = SS_update(ss_table, diff_t, buses, dct_flag, bus_e, t_updt, index)
            t = t_updt
            prev_time[index] = t
            
            # if its the last event for the bus, update bus parameters to standstill
            if buses[index].next_e()==np.inf:
                buses[index].state = -1
                buses[index].route = None
                buses[index].time_arr = list()
                buses[index].event_arr = list()
                prev_time[index] = np.nan
            index = np.nan
            t_updt = np.nan
            bus_e = np.nan

        #-----
        # Case 7: Jumping to next demand 
        elif ((t>=T) and (bus_chk==-1) and
              (demand_ct[0]==time_check) and 
              (demand_ct[0]!=np.inf)):
            #print('\tJumping to next demand')
            msg = 'Jumping to next demand'
            t = next_bus_e(buses)[0]
            #dct_flag += 1
        
        #-----
        # Case 8: All demands are completed outside the timeframe 
        elif ((t>=T) and
              (time_check==np.inf)):

            #print('\tJumping to EOD')
            msg = 'Jumping to EOD'
            t_updt = np.inf        
            ss_table = SS_update(ss_table, t, buses, dct_flag, bus_e, t_updt)
            t = t_updt    
        
        # updating all current demand times
        demand_ct = list(np.sort(demand_ct))
        for k in range(dct_flag):
            demand_ct[k] = t
            
    return ss_table, bd_table


def cost_analysis(n_buses, replicates, refuel, ss_table, bd_table, emp_rate, fuel_rate, delay_rate, refuel_consumption):
    ss_table['Process_Time'].replace({np.inf: 0}, inplace=True)
    ss_table['Demand_Charge'].replace({np.nan: -1}, inplace=True)
    cost_ss = ss_table[(ss_table['Demand_Charge'])==-1][['Bus', 'Route', 'Event', 'Process_Time']].dropna()

    runcost_ss = cost_ss[cost_ss['Event']==1]
    sercost_ss = cost_ss[(cost_ss['Event']!=1)&(cost_ss['Route']!=refuel)]
    refcost_ss = cost_ss[(cost_ss['Event']!=1)&(cost_ss['Route']==refuel)]
    runcost_ss = runcost_ss[['Bus', 'Process_Time']].groupby('Bus').sum()/replicates
    sercost_ss = sercost_ss[['Bus', 'Process_Time']].groupby('Bus').sum()/replicates
    refcost_ss = refcost_ss[['Bus', 'Process_Time']].groupby('Bus').sum()/replicates

    run_time = round(np.sum(runcost_ss['Process_Time']), 2)
    ser_time = round(np.sum(sercost_ss['Process_Time']), 2)
    ref_time = round(np.sum(refcost_ss['Process_Time']), 2)
    emp_cost = round((run_time + ser_time + ref_time)*emp_rate/60)
    fuel_cost = round(ref_time*refuel_consumption*-1*fuel_rate)
    tot_delay = round(sum(np.nan_to_num(np.array(bd_table['Demand_Current'] - bd_table['Demand_Actual']))), 2)
    delay_cost = round(tot_delay* delay_rate)
    return run_time, ser_time, ref_time, emp_cost, fuel_cost, delay_cost


# +
def ecdf(target, title):
    numbs = np.array(target)
    ecdf = sm.distributions.ECDF(numbs)
    x = np.linspace(min(numbs), max(numbs), len(target))
    y = ecdf(x)
    plt.step(x, y, color='r')
    plt.xlabel(title)
    plt.ylabel('CDF(P)')
    plt.title('Empirical CDF of ' + title)
    plt.show()
    return

def inv_ecdf(target, title):
    numbs = np.array(target)
    ecdf = sm.distributions.ECDF(numbs)
    x = np.linspace(min(numbs), max(numbs), len(target))
    y = 1-ecdf(x)
    plt.step(x, y, color='r')
    plt.xlabel(title)
    plt.ylabel('Inverse CDF(P)')
    plt.title('Inverse Empirical CDF of ' + title)
    #plt.gca().invert_yaxis()
    plt.show()
    return
