#!/home/anaconda/bin/python
import pandas as pa
import pylab
import numpy as np
from numpy import arange, argwhere, median
import argparse


if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='find where derivative of data os ')
	parser.add_argument('-i','--input',default='', help='input file')
	parser.add_argument('-c','--conservative',  action='store_true' , help='conservative')
	parser.add_argument('-p','--plot',  action='store_true' , help='plot')
	args = parser.parse_args()

	#y =  np.genfromtxt(args.input)
	t=pa.read_csv (args.input,sep='\t',header=None)
	#y=t[3] - t[4]
	y=t[7]
	y-=min(y)
	y/=max(y)

	xdata = arange(len(y))*1. / len(y)
	smoothed_y  = np.array(pa.ewma(y,span=100)).ravel()
	diff_y = np.array(pa.ewma(np.diff(smoothed_y),span=100)).ravel() / (xdata[1]-xdata[0])
#	print min(argwhere(np.diff(smoothed_y)*len(y)>1).ravel())
	if args.conservative :
		print argwhere(np.diff(diff_y >1)).ravel()[-1]
	else:
		print int(median(argwhere(np.diff(diff_y >1)).ravel() ))
	if args.plot:
		pylab.scatter(xdata,y,marker='.',s = 1)
		pylab.plot(xdata,smoothed_y)
		pylab.scatter(xdata[1:len(y)] , np.diff(smoothed_y)/np.diff(xdata)[0]>1 ,s=.5)
		pylab.savefig('out.pdf')	
