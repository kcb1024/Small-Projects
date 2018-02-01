#!/usr/bin/python3
import sys
import os
import argparse
import pandas
import matplotlib.pyplot as plt

# Get the version of python being used to run the script
if (sys.version_info[0]!=3 ):
    print(sys.version_info)
    raise Exception("Run with python3")

# Parse arguments
parser = argparse.ArgumentParser()
parser.add_argument("-f", help="specify file to read", type = str)

args = parser.parse_args()
filepath = args.f
print(filepath)

# Read in the data with pandas
df = pandas.read_table(str(filepath), sep=' |=|:', header = None, engine = 'python')

# Check to see if Humidity, Temp strings added to array
dotw = ['Sun', 'Mon', 'Tues', 'Wed', 'Thurs', 'Fri', 'Sat']
undesirables = ['Temp', 'Humidity']
#undesirables += dotw

# Screen for unwanted data by examining first row of data
dat = df[:][0:1]
dat_len = dat.size
undesirable_idx = []
for u in undesirables:
    for i in range(dat_len):
        a = df.at[0,i]
        # If one of the undesired fields in dataframe, add column idx to list
        if (a == u):
            undesirable_idx += [i]

# Remove unwanted data
df = df.drop(columns = undesirable_idx)

# Rename column names
num_cols = df[:][0:1].size
# Get current column names
c = df.columns

# This is the one implemented already
if (num_cols == 10):
    cols = ['Temp', 'Humidity', 'DOTW', 'Day', 'Month', 'Hour', 'Minutes', 'Seconds', 'Timezone', 'Year']
    col_dict = dict()
    for old_col in c:
        col_dict[old_col] = cols[0]
        # Remove zeroth entry, since it's already been used
        cols = cols[1:]
    df = df.rename(index = str, columns = col_dict)

# Clean up unwanted characters in temp (*) and humidity (%) fields and convert to floats
df['Temp'] = df['Temp'].map(lambda x: x.rstrip('*') )
df['Temp'] = df['Temp'].astype(float)
df['Humidity'] = df['Humidity'].map(lambda x: x.rstrip('%') )
df['Humidity'] = df['Humidity'].astype(float)

print(df.head() )
temp = df['Temp'].as_matrix()
humidity = df['Humidity'].as_matrix()
plt.figure()
plt.plot(temp, humidity, 'b.')
plt.show()
