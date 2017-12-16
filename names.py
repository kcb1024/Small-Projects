# -*- coding: utf-8 -*-
"""
Created on Fri Dec 08 22:05:40 2017

@author: karen
"""

import numpy as np
import os
import time
import matplotlib.pyplot as plt
import random


"""
Designed to estimate the age distribution of a group based off the names
of the people in the group. It uses data from the Social Security
Administration to estimate age based off the frequency that first names
occur in the data, because the popularity of names changes over time.

Example commands:
# Creates object scanning SSA data located in dir_path over specified time range
n = Names(dir_path),1950,2000, 0.01)
# Add in a list of names (>100 names is best) for analysis
l = ['John', 'Anna', 'Henry', 'Ophelia', 'Henry', 'Allison', 'Allison']
# Analyze the list
n.create_prob_mapping(l)
"""


class Names:
    def __init__(self, 
                 dir_path,
                 start_year,
                 stop_year,
                 inclusion_threshold_pct = 0.01):
        self.dir_path = dir_path
        self.start_year = start_year
        self.stop_year = stop_year
        self.name_idx = {}
        self.num_total_births = []
        self.inclusion_threshold_pct = inclusion_threshold_pct

        # Row = name index
        # Col = year
        self.name_data = []
        
        # Run data
        self.read_ssa_data()
        
    def read_ssa_data(self):
        """
        Gets data from set of Social Security Administration 
        files and stores in self.name_idx and self.name_data.
        """
        t0 = time.time()
        
        # First, build list of names to be evaluated

        # Get list of all files from SSA
        files_list_temp = os.listdir(self.dir_path)
        # Trim list of files to only desired years
        files_list = []
        # Strip year out of file name, e.g. 'yob2016.txt'
        s1 = 3
        s2 = 7
        # Add any files inside the desired range to the new list of files
        for f in files_list_temp:
            if f[0:3] == 'yob':
                if self.start_year <= int(f[s1:s2]) <= self.stop_year:
                    files_list += [f]
        # Make sure dir path ends with slash
        if self.dir_path[-1] != '/':
            self.dir_path += '/'
        # Scan through files for names to add
        ssa_idx_name = 0
        ssa_idx_num_births = 2
        ctr = 0        
        
        for filename in files_list:
            total_annual_births = 0
            # Read in data
            f = open(self.dir_path + filename)
            data = f.read()
            ds = data.splitlines()
            # Calculate total number of births
            for i in ds:
                i = i.split(',')
                total_annual_births += int(i[ssa_idx_num_births])
            self.num_total_births += [total_annual_births]

            # Go through and add popular names
            # Calculate pct for all names, changing raw number -> pct in array
            for i in ds:
                i = i.split(',')
                pct = int(i[ssa_idx_num_births])/float(total_annual_births)*100
                if pct*2 > self.inclusion_threshold_pct:
                    # check to see if the name is already listed
                    if i[ssa_idx_name] not in self.name_idx:
                        # add it
                        name = i[ssa_idx_name]
                        self.name_idx[name] = ctr
                        # increment counter
                        ctr += 1
#        print(total_num_births_list)
        # build array of popular names over time
        num_rows = len(self.name_idx)
        num_cols = self.stop_year - self.start_year + 1
        self.name_data = np.zeros([num_rows, num_cols])
        yr_ctr = 0
        for filename in files_list:
            # Read in data
            f = open(self.dir_path + filename)
            data = f.read()
            ds = data.splitlines()
            # Go through each name
            for j in ds:
                # split one long string into multiple strings
                i = j.split(',')
                # Check to see if the name is popular enough to care about it
                if i[ssa_idx_name] in self.name_idx:
                    name = i[ssa_idx_name]
                    idx = self.name_idx[name]
                    pct = float(i[ssa_idx_num_births])/self.num_total_births[yr_ctr]*100
                    self.name_data[idx,yr_ctr] += pct                
            yr_ctr+=1
        t1 = time.time()
        print('time elapsed: ', t1 - t0)
        return
        
    def build_prob_map_one_name(self, name, name_prop, lower_lim):
        """ 
        Given the proportion (pct) of a name in population, generate map of
        probabilities that group was born in given year
        """
        yrs = range(self.start_year, self.stop_year+1)
        mapping = np.zeros(len(yrs))
        # check that name is in the registry
        if name not in self.name_idx:
            return []
        # scan through years
        for y in yrs:
            y_i = y - self.start_year
            theo = self.name_data[self.name_idx[name],y_i]
            if (theo > 0):
                experimental = name_prop
                error = abs(theo - experimental)/theo
    #            mapping[y_i] = (1 - (error)**2)*theo #1            
    #            mapping[y_i] = (1 - (error)**2)*theo #2
                mapping[y_i] = (error*theo)**2 #3
    #            mapping[y_i] *= (lower_lim < mapping[y_i] < 1)
            
        return mapping
    def create_prob_mapping(self, names_list, lower_lim = 0):
        """
        Takes in a list of names.
        Spits out a plot of overlayed probability maps for each name
        """
        plt.figure()
        plt.subplot(2,1,1)
        # convert from raw list to a dict
        sample_pop = convert_simple_list_to_dict(names_list)
        # convert from raw numbers to proportions
        sample_pop = convert_names_to_props(sample_pop)
#        print(sample_pop)
        yrs = range(self.start_year, self.stop_year + 1)
        all_names = np.zeros(len(yrs))
        for name in sample_pop.keys():
            prop = sample_pop[name]
            mapping = self.build_prob_map_one_name(name, prop, lower_lim)
#            print(mapping)
            if mapping != []:
                all_names += mapping
                plt.plot(yrs, mapping, label = name)
        plt.xlabel('years')
        plt.subplot(2,1,2)
        plt.plot(yrs, all_names, 'b.-')
        return yrs, all_names
#        plt.legend()
    def calc_pct_popular_names(self, names_list):
        """
        Input: 1-d list of names
        Output: pct of entries in list that are in popular names list        
        """
        num_names_popular = 0
        for name in names_list:
            if name in self.name_idx:
                num_names_popular += 1
        pct = float(num_names_popular)/len(names_list)*100
        return pct
    def build_sample_population(self, yr_dict, pop_size):
        """
        Generates a sample population based off SSA data.
        Useful for verification of func create_prob_mapping
        Input:
        - yr_dist: dict with pct of population born in each year
        - pop_size (int): total number of people in desired sample population
        Output:
        - sample population: 1-d list of names of people in population
        """
        # create empty list for output population
        sample_pop = []
        # iterate through each year that sample pop is born in
        # first, normalize the values to 100%, b/c user might have been lazy
        # and not made the pct add up to 100%
        total_pct = 0
        for yr in yr_dict.keys():
            total_pct += yr_dict[yr]
        # now go forth and start generating names
        for yr in yr_dict.keys():
#            print(yr)
            # normalize pct to 100%
            pct_yr = yr_dict[yr]/total_pct*100
            # scan through SSA data and assign idx to names to reflect
            # popularity, so we can randomly pick them according to their
            # popularity
            # The problem with this is that is assumes popular names are 100%
            # of the names in the data set
            yr_pop_arr = []
            pop_ctr = 0
            for name in self.name_idx:
                pop_ctr += self.name_data[self.name_idx[name], yr - self.start_year]
                # 0 column is the population fraction (0 - 1)
                yr_pop_arr += [pop_ctr]
                # increment counter for the array index
#                idx_ctr += 1
            # scale by number of total births that year
            yr_pop_arr = np.array(yr_pop_arr)*1.0
            yr_pop_arr /= yr_pop_arr[-1]
            # now start generating names
            # get number of names for this year
#            print('total births this year ', yr_pop_arr[-1])
            num_names_yr = int(round(pct_yr/100.0*pop_size))
#            print(num_names_yr)
            for i in range(num_names_yr):
                # generate random number in [0,1]
                r = random.random()
                idx = sort_to_closest(yr_pop_arr, r)
                name_items = self.name_idx.items()
                name = name_items[idx][0]
                sample_pop += [name]
                print(r, idx, name, self.name_data[self.name_idx[name], yr - self.start_year]/yr_pop_arr[-1])
        
        return sample_pop

def sort_to_closest(l, r):
    """
    takes in a sorted list and returns the index of the closest
    value to r
    """
    # find closest value to r
#    ltr = 1 # less than r
    idx = -1
    ctr = 0
    for i in l:
        # if i not less than r
        if (i > r):
            return ctr
        ctr += 1
    return idx

        
def convert_names_to_props(sample_pop):
    """ takes in a dict 
    spits out proportions in pct    
    """
    total_pop = 0
    # calculate total population
    keys = sample_pop.keys()
    for name in keys:
        total_pop += float(sample_pop[name])
    # calculate individual proportions
    for name in keys:
        sample_pop[name] = (sample_pop[name]/total_pop)*100
    return sample_pop
def convert_simple_list_to_dict(l):
    """
    takes 1-d list of names
    generates a dict
    """
    d = {}
    for name in l:
        if name not in d:
            d[name] = 1
        else:
            d[name] += 1
    return d