##################################
# Fiducial intra rater calculator #
##################################

import sys
import numpy as np

# importing csv module
import csv

# csv file name
filename = sys.argv[1]
 
# initializing the titles and rows list
fields = []
rows = []
 
# reading csv file
with open(filename, 'r') as csvfile:
    # creating a csv reader object
    csvreader = csv.reader(csvfile)

    # extracting field names through first row
    comment1 = csvreader.next()
    comment2 = csvreader.next()
    fields = csvreader.next()
 
    # extracting each data row one by one
    for row in csvreader:
        rows.append(row)

x_values = []
y_values = []
z_values = []

for current_row in rows:
	x_values.append(current_row[1])
	y_values.append(current_row[2])
	z_values.append(current_row[3])

# print x_values
# print y_values
# print z_values