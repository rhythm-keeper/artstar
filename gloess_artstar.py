#!/usr/bin/env

import os
import numpy as np
from astropy.io import fits
import random
import bisect


def file_len(fname):
        with open(fname) as f:
                for i, l in enumerate(f):
                        pass
        return i + 1


def make_ascii_record(artstar_file):

	# Read data into columns first,
	# then store in a records table for cleanness 

	n_stars_max=file_len(artstar_file)-1 # -1 to ignore heaer

	record = np.zeros( n_stars_max, dtype=[('x_in','f8'),('y_in','f8'),('f606w_in','f8'),\
						('f814w_in','f8'),('f606w_out','f8'),('f814w_out','f8')] )

	rownum=0
	with open(artstar_file) as f:

		header=f.readline()

		for row in f:

			parts=row.split()

			record[rownum]=(float(parts[1]),float(parts[2]),float(parts[3]),float(parts[4]),\
				float(parts[5]),float(parts[7]))

			rownum+=1

	return record


def get_artstar_subsample(record_table,nstars_to_sample,use_cdfs=0):


	num_poss=len(record_table)
	

	# Ext 1 chip 2
	subset_indices=np.zeros(nstars_to_sample, dtype=int)
	#x_artstar=np.zeros(nstars_to_sample)
	#y_artstar=np.zeros(nstars_to_sample)
	#f814w_artstar=np.zeros(nstars_to_sample)
	#f606w_artstar=np.zeros(nstars_to_sample)
                                             

	count=0
	while count < nstars_to_sample: 
	                                                                                                                             
		if use_cdfs: 
	                                                                                                                             
			rand_x_val=randomu(seed)
			rand_y_val=randomu(seed)
			
			# find where the random number falls in the CDF
			prob_x_cdf_val_of_nearest_result=nearest_element(rand_x_val,prob_x_cdf,ind_of_nearest_x_el)
			#dx=abs(prob_x_cdf-rand_x_val)
			#ind_of_nearest_x_el=where(dx eq min(dx))
			
			prob_y_cdf_val_of_nearest_result=nearest_element(rand_y_val,prob_y_cdf,ind_of_nearest_y_el)	
			# dy=abs(prob_y_cdf-rand_y_val)
			# ind_of_nearest_y_el=where(dy eq min(dy))
			
			# find out what x value this place in the CDF corresponds to
			xval_of_nearest_el=xpos[ind_of_nearest_x_el]              			
			yval_of_nearest_el=ypos[ind_of_nearest_y_el]
			
			# add random element to positions to pick out 
			# so the same sources don't get reselect every time
			xval_of_nearest_el+=randomn(seed)*rwalk_scale # random walk	
			yval_of_nearest_el+=randomn(seed)*rwalk_scale #
			
			
			search_rad = 10
			match = -1
			while match == -1: 
			        match=match_2d([xval_of_nearest_el],[yval_of_nearest_el],x_in,y_in,search_rad,MATCH_DISTANCE=md)
			        #print,xval_of_nearest_el,yval_of_nearest_el,x_in[match],y_in[match]
			        search_rad *= 2
			#print,xval_of_nearest_el,x_in[match],yval_of_nearest_el,y_in[match]
			
			#el=0L
			#while el lt n_elements(x_in) do begin
			#	diff_x_vals[el]=x_in[el]-xval_of_nearest_el
			#	diff_y_vals[el]=y_in[el]-yval_of_nearest_el
			#	el+=1
			#endwhile 
				
			# find closest xval in simulated data
			 #dx=x_in-xval_of_nearest_el
			 #dy=y_in-yval_of_nearest_el
			 #distances=sqrt(diff_x_vals^2.+diff_y_vals^2.)
			 #min_dist_ind=where(distances eq min(distances))
			
			
			
			#print,strtrim(xval_of_nearest_el,2)+" "+strtrim(x_in[min_dist_ind],2)+" "+$
			#	strtrim(yval_of_nearest_el,2)+" "+strtrim(y_in[min_dist_ind],2)
		else:
			# random selection
			rand_int=random.randint(0,num_poss-1)
			subset_indices[count]=rand_int	                                                                                                                             

	
		# keep values of artstar

		#x_artstar[count]=record_table['x_in'][match]
		#y_artstar[count]=record_table['y_in'][match]
		#f814w_artstar[count]=record_table['f814w_out'][match]
		#f606w_artstar[count]=record_table['f606w_out'][match]
		 #agb check
		#if id_artstar[match] ge 3e5 then is_agb[count]=1 else is_agb[count]=0
	                                                                                                                             
		count+=1


	return subset_indices#[f606w_artstar, f814w_artstar]


def count_num_in_bin(data,bin_el,binsize,min_data_val):

	# Get starting point to count points
	start_loc=bisect.bisect_left(data,min_data_val+binsize*bin_el) 

	low_bound   = min_data_val + binsize * bin_el
	upper_bound = min_data_val + binsize * (bin_el + 1)

	# Find first valid star in range
	# Possible bisect returns first
	# valid star immediately if
	# data val == low_bound
	while data[start_loc] < low_bound:
		start_loc += 1


	# Now count the rest in the range
	num_star=0
	el = start_loc
	max_el = len(data)-1
	while data[el] <= upper_bound:
		num_star += 1
		el += 1
		if el >= max_el: break # prevent indexing error

	return num_star



def gloess(data,binsize,smooth_factor):

	#==========================
	# First bin the data based
	# on the smooth_factor
	#==========================

	# sort data to for faster processing
	ind=data.argsort()
	data=data[ind]

	# Set boundaries for LF
	max_data_val=np.max(data)
	min_data_val=np.min(data)

	n_data_bins=int(np.ceil( (max_data_val-min_data_val)/binsize ))
	
	binned_array_count=np.zeros(n_data_bins)
	binned_array_value=np.zeros(n_data_bins)

	for bin_el in range(n_data_bins):

		# Get number within mag bin range	
		num_in_bin = count_num_in_bin(data,bin_el,binsize,min_data_val)
	
		binned_array_count[bin_el] = num_in_bin
		binned_array_value[bin_el] = min_data_val + binsize * (bin_el + 0.5)
	
	
	#=============================================
	# Now loop through an evenly spaced
	# grid of points based on the smooth_factor, 
	# compute the distance between points, compute
	# the weighted contribution based on a Gaussian
	# then record the weighted mean of all points
	#=============================================
	
	weight=np.zeros(n_data_bins) # for containing a Gaussian weights from center of bin
	new_LF=np.zeros(n_data_bins) # to contained smoothed bin

	for bin_el in range(n_data_bins):

		# The position used as the current reference point
		curr_position = min_data_val + binsize * (bin_el + 0.5)
	
		# Loop through other data points and compute their contribution
		for other_bin_el in range(n_data_bins):
	
			other_position = min_data_val + binsize * (other_bin_el + 0.5)

			distance_from_curr_pos = abs(curr_position - other_position)
	        
			exp_fact = -0.5 * (distance_from_curr_pos / smooth_factor)**2.
			weight[other_bin_el] = np.exp( exp_fact )
	
	
	        # Compute weighted mean
		dum=0
		total_weight = np.sum(weight)
		for bin_loop in range(n_data_bins):
			frac_weight = weight[bin_loop] / total_weight
			dum += binned_array_count[bin_loop] * frac_weight 

		new_LF[bin_el]=dum
	
	
	
	# Compute Sobel [-1,0,+1] edge response and S/N function
	sobel   = np.zeros(n_data_bins-2)
	sigma   = np.zeros(n_data_bins-2)
	bin_val = np.zeros(n_data_bins-2)
	rownum = 0
	for bin_loop in range(1,n_data_bins-1):
		sobel[rownum] = new_LF[bin_loop+1] - new_LF[bin_loop-1] 
		sigma[rownum] = np.sqrt( new_LF[bin_loop+1] + new_LF[bin_loop-1] ) 
		bin_val[rownum] = min_data_val + binsize * (bin_loop + 0.5)
		rownum += 1
	
	total_weight = np.sum(weight)
	frac_weight = [abs(x) / y * 1./total_weight for x,y in zip(sobel,sigma)]
	                                                                                               
	weighted_sobel = [x * y for x,y in zip(sobel,frac_weight)]
	
	# get location of maximum in edge response
	# WHAT TO DO IF EQUAL PEAKS?
	ind_of_max = weighted_sobel.index(np.max(weighted_sobel))	

	return bin_val[ind_of_max]





def main():


	# spatially sample artificial stars
	# to follow distribution in images
	use_cdfs=0
	# if using cdfs, set random walk scale
	# to use other nearby stars at random
	# if the sample of artificial stars
	# is small. other same stars may be
	# frequently chosen
	rwalk_scale=10 # pix
	
	# the range in mag centered
	# on the trgb to sample stars from
	mag_range=0.75
	
	# number of trgb measurements
	n_attempts=int(1e4)
	
	# if the dao output
	# hasnt been calibrated yet
	# put it on a scale to make
	# the trgb simulation output
	# seem familiar
	F814W_cal=7.2
	input_TRGB=26.3
	
	
	# read in artificial stars data for both chips
	record_table_1 = make_ascii_record("condensed_artstar_1.dat") 
	record_table_4 = make_ascii_record("condensed_artstar_4.dat")

	#print(record_table_4['x_in'][:])
	#exit(0)
	
	
	#if use_cdfs:
	#	# CDFs to sample stars
	#	# Convert simulation uniform X-Y to prob distr
	#	n_stars_max=file_len("NGC1316_x_cdf_1.dat")
	#	xpos_1=np.zeros(n_stars_max)
	#	prob_x_cdf_1=np.zeros(n_stars_max)
	#	n_stars_max=file_len("NGC1316_x_cdf_1.dat")
	#	xpos_4=np.zeros(n_stars_max)
	#	prob_x_cdf_4=np.zeros(n_stars_max)
	#	n_stars_max=file_len("NGC1316_y_cdf_1.dat")
	#	ypos_1=np.zeros(n_stars_max)
	#	prob_y_cdf_1=np.zeros(n_stars_max)
	#	n_stars_max=file_len("NGC1316_y_cdf_4.dat")
	#	ypos_4=np.zeros(n_stars_max)
	#	prob_y_cdf_4=np.zeros(n_stars_max)
	
	
	
	
	min_smoothing_scale=0.09
	max_smoothing_scale=0.20
	dsmoothing_scale=0.01
	binsize=0.01
	nsmoothing_scales=int((max_smoothing_scale-min_smoothing_scale)/dsmoothing_scale)
	for smoothing_scale_loop in range(nsmoothing_scales):
	
		
		

		# Get catalog of artificial stars
		# read catalog
		
		smooth_factor = smoothing_scale_loop * dsmoothing_scale + min_smoothing_scale
	
		
		# Loop through many GLOESS fit attempts
		
		#output_label=str(smooth_factor,format='(F5.3)')
		output_label=str(("{0:.3f}").format(smooth_factor))
		os.system("rm -f ./artstar_trgb_output/artstar_gloess_edge_detections_"+output_label+".dat")
		with open("./artstar_trgb_output/artstar_gloess_edge_detections_"+output_label+".dat","a") as outf:


			
			

		
			# reduce by factor is the number within the color-mag cut compared to the entire sample
			# in the mag range of the trgb
			nstars_to_sample_1 = 607 # +/-0.75 mag
			nstars_to_sample_4 = 166 # +/-0.75 mag	
	
			for gloess_attempt in range(n_attempts):
	
				# Get stars chip 1
				subset_indices_4 = get_artstar_subsample(record_table_4, nstars_to_sample_4)
				# Get stars chip 2
				subset_indices_1 = get_artstar_subsample(record_table_1, nstars_to_sample_1)
	
				# Combine stars from both CCDs
				f814w_artstar=np.array( list(record_table_1['f814w_out'][subset_indices_1]) +\
							list(record_table_4['f814w_out'][subset_indices_4]) ) 

				measured_trgb = gloess(f814w_artstar,binsize,smooth_factor)

				print(measured_trgb)
	
				outf.write(str(measured_trgb)+"\n")

		


if __name__ == "__main__": main()
