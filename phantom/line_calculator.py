##################################################
# Slicer Fiducial Point Line Distance Calculator #
##################################################

import sys
import numpy as np

# importing csv module
import csv


# csv file name
filename = sys.argv[1]

#number of fiducials
 
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


## Left ##

pt = np.array([float(x_values[7]), float(y_values[7]), float(z_values[7])])
x1 = np.array([float(x_values[0]), float(y_values[0]), float(z_values[0])])
x2 = np.array([float(x_values[1]), float(y_values[1]), float(z_values[1])])

# pt = np.array([1.0, 0, 1.0])
# x1 = np.array([1.0, 2.0, -1.0])
# x2 = np.array([2.0, 0, 3.0])

x1_minus_pt = pt - x1

x2_minus_x1 = x2 - x1

sumsq_x1_minus_pt = sum(x1_minus_pt * x1_minus_pt)

sumsq_x2_minus_x1 = sum(x2_minus_x1 * x2_minus_x1)

mydotprod = np.dot(x1_minus_pt, x2_minus_x1)

dist3d = np.sqrt((sumsq_x1_minus_pt * sumsq_x2_minus_x1 - (mydotprod * mydotprod))/sumsq_x2_minus_x1)

# print "pt: ", pt
# print "x1: ", x1
# print "x2: ", x2
# print "x1_minus_pt: ", x1_minus_pt
# print "x2_minus_x1: ", x2_minus_x1
# print "sumsq_x1_minus_pt: ", sumsq_x1_minus_pt
# print "sumsq_x2_minus_x1: ", sumsq_x2_minus_x1
# print "mydotprod: ", mydotprod
print "Left : ", dist3d

## Right ##

pt = np.array([float(x_values[8]), float(y_values[8]), float(z_values[8])])
x1 = np.array([float(x_values[0]), float(y_values[0]), float(z_values[0])])
x2 = np.array([float(x_values[2]), float(y_values[2]), float(z_values[2])])

# pt = np.array([1.0, 0, 1.0])
# x1 = np.array([1.0, 2.0, -1.0])
# x2 = np.array([2.0, 0, 3.0])

x1_minus_pt = pt - x1

x2_minus_x1 = x2 - x1

sumsq_x1_minus_pt = sum(x1_minus_pt * x1_minus_pt)

sumsq_x2_minus_x1 = sum(x2_minus_x1 * x2_minus_x1)

mydotprod = np.dot(x1_minus_pt, x2_minus_x1)

dist3d = np.sqrt((sumsq_x1_minus_pt * sumsq_x2_minus_x1 - (mydotprod * mydotprod))/sumsq_x2_minus_x1)

# print "pt: ", pt
# print "x1: ", x1
# print "x2: ", x2
# print "x1_minus_pt: ", x1_minus_pt
# print "x2_minus_x1: ", x2_minus_x1
# print "sumsq_x1_minus_pt: ", sumsq_x1_minus_pt
# print "sumsq_x2_minus_x1: ", sumsq_x2_minus_x1
# print "mydotprod: ", mydotprod
print "Right : ", dist3d

## Superior ##

pt = np.array([float(x_values[9]), float(y_values[9]), float(z_values[9])])
x1 = np.array([float(x_values[0]), float(y_values[0]), float(z_values[0])])
x2 = np.array([float(x_values[3]), float(y_values[3]), float(z_values[3])])

x1_minus_pt = pt - x1

x2_minus_x1 = x2 - x1

sumsq_x1_minus_pt = sum(x1_minus_pt * x1_minus_pt)

sumsq_x2_minus_x1 = sum(x2_minus_x1 * x2_minus_x1)

mydotprod = np.dot(x1_minus_pt, x2_minus_x1)

dist3d = np.sqrt((sumsq_x1_minus_pt * sumsq_x2_minus_x1 - (mydotprod * mydotprod))/sumsq_x2_minus_x1)

# print "pt: ", pt
# print "x1: ", x1
# print "x2: ", x2
# print "x1_minus_pt: ", x1_minus_pt
# print "x2_minus_x1: ", x2_minus_x1
# print "sumsq_x1_minus_pt: ", sumsq_x1_minus_pt
# print "sumsq_x2_minus_x1: ", sumsq_x2_minus_x1
# print "mydotprod: ", mydotprod
print "Superior : ", dist3d

## Inferior ##

pt = np.array([float(x_values[10]), float(y_values[10]), float(z_values[10])])
x1 = np.array([float(x_values[0]), float(y_values[0]), float(z_values[0])])
x2 = np.array([float(x_values[4]), float(y_values[4]), float(z_values[4])])

x1_minus_pt = pt - x1

x2_minus_x1 = x2 - x1

sumsq_x1_minus_pt = sum(x1_minus_pt * x1_minus_pt)

sumsq_x2_minus_x1 = sum(x2_minus_x1 * x2_minus_x1)

mydotprod = np.dot(x1_minus_pt, x2_minus_x1)

dist3d = np.sqrt((sumsq_x1_minus_pt * sumsq_x2_minus_x1 - (mydotprod * mydotprod))/sumsq_x2_minus_x1)

# print "pt: ", pt
# print "x1: ", x1
# print "x2: ", x2
# print "x1_minus_pt: ", x1_minus_pt
# print "x2_minus_x1: ", x2_minus_x1
# print "sumsq_x1_minus_pt: ", sumsq_x1_minus_pt
# print "sumsq_x2_minus_x1: ", sumsq_x2_minus_x1
# print "mydotprod: ", mydotprod
print "Inferior : ", dist3d

## Anterior ##

pt = np.array([float(x_values[11]), float(y_values[11]), float(z_values[11])])
x1 = np.array([float(x_values[0]), float(y_values[0]), float(z_values[0])])
x2 = np.array([float(x_values[5]), float(y_values[5]), float(z_values[5])])

x1_minus_pt = pt - x1

x2_minus_x1 = x2 - x1

sumsq_x1_minus_pt = sum(x1_minus_pt * x1_minus_pt)

sumsq_x2_minus_x1 = sum(x2_minus_x1 * x2_minus_x1)

mydotprod = np.dot(x1_minus_pt, x2_minus_x1)

dist3d = np.sqrt((sumsq_x1_minus_pt * sumsq_x2_minus_x1 - (mydotprod * mydotprod))/sumsq_x2_minus_x1)

# print "pt: ", pt
# print "x1: ", x1
# print "x2: ", x2
# print "x1_minus_pt: ", x1_minus_pt
# print "x2_minus_x1: ", x2_minus_x1
# print "sumsq_x1_minus_pt: ", sumsq_x1_minus_pt
# print "sumsq_x2_minus_x1: ", sumsq_x2_minus_x1
# print "mydotprod: ", mydotprod
print "Anterior : ", dist3d

## Posterior ##

pt = np.array([float(x_values[12]), float(y_values[12]), float(z_values[12])])
x1 = np.array([float(x_values[0]), float(y_values[0]), float(z_values[0])])
x2 = np.array([float(x_values[6]), float(y_values[6]), float(z_values[6])])

x1_minus_pt = pt - x1

x2_minus_x1 = x2 - x1

sumsq_x1_minus_pt = sum(x1_minus_pt * x1_minus_pt)

sumsq_x2_minus_x1 = sum(x2_minus_x1 * x2_minus_x1)

mydotprod = np.dot(x1_minus_pt, x2_minus_x1)

dist3d = np.sqrt((sumsq_x1_minus_pt * sumsq_x2_minus_x1 - (mydotprod * mydotprod))/sumsq_x2_minus_x1)

# print "pt: ", pt
# print "x1: ", x1
# print "x2: ", x2
# print "x1_minus_pt: ", x1_minus_pt
# print "x2_minus_x1: ", x2_minus_x1
# print "sumsq_x1_minus_pt: ", sumsq_x1_minus_pt
# print "sumsq_x2_minus_x1: ", sumsq_x2_minus_x1
# print "mydotprod: ", mydotprod
print "Posterior : ", dist3d
