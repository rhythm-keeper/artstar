#!/usr/bin/env


import os
import numpy as np
import glob
import astropy.io.fits as pyfits

import bisect


def file_len(fname):
	with open(fname) as f:
		for i, l in enumerate(f):
			pass
	return i + 1

def get_avg_mag_magerr(mag,magerr):
	"""
	Returns flux averaged magnitudes
	Inputs: mag array, magerr array
	"""
	
	fluxes=[]
	fluxerrs=[]
	
	for mag_el,magerr_el in zip(mag,magerr):
	        flux_el=10**( (mag_el-25.)/(-2.5) )
	        fluxerr_el=flux_el*magerr_el/1.0857
	
	        fluxes.append(flux_el)
	        fluxerrs.append(fluxerr_el)
	
	# wavg flux
	tot_weight=0
	dum=0
	for flux_el,fluxerr_el in zip(fluxes,fluxerrs):
	        weight=1./fluxerr_el**2
	        tot_weight+=1./fluxerr_el**2
	        dum+=flux_el*weight
	
	avg_flux=dum/tot_weight
	
	inv_sqr_flux=[1./x**2 for x in fluxerrs]
	avg_fluxerr=np.sqrt( 1./np.sum(inv_sqr_flux) )
	
	avg_mag=-2.5*np.log10(avg_flux)+25
	avg_magerr=avg_fluxerr*1.0857/avg_flux
	
	return [avg_mag,avg_magerr]


def get_artstar_coord_to_match_against_raw_file(f814w_add_file):

	nrows=file_len(f814w_add_file)
	x_artstar=np.zeros(nrows-3)
	y_artstar=np.zeros(nrows-3)
	mag_artstar=np.zeros(nrows-3)
	rownum=0
	with open(f814w_add_file,"r") as f:

		# Ignore header
		for x in range(3): f.readline()

		for row in f:
			parts=row.split()
			x_artstar[rownum]=parts[1]
			y_artstar[rownum]=parts[2]
			mag_artstar[rownum]=parts[3]
	                                                                                                  
			rownum+=1
	                                                                                                  
	# Sort by x_array
	inds=x_artstar.argsort()
	x_artstar=x_artstar[inds]
	y_artstar=y_artstar[inds]
	mag_artstar=mag_artstar[inds]
	#print(x_artstar[0:3],y_artstar[0:3])
	


	return [x_artstar,y_artstar,mag_artstar]


def load_pull_write_artstars(raw_file,f814w_add_file,f606w_add_file,filebasenames,exptime_mag_diff,el,n_images,ref_f814w_exptime,ref_f606w_exptime,n_f814w_images):


	# get coordinate list of input artificial stars
	mag1_x_artstar,mag1_y_artstar,mag1_artstar=get_artstar_coord_to_match_against_raw_file(f814w_add_file)
	mag2_x_artstar,mag2_y_artstar,mag2_artstar=get_artstar_coord_to_match_against_raw_file(f606w_add_file)

	
	
	# count number of stars in file
	num_rows_per_object=int(n_images / 6) + 1 # 6 images per line
	nlines=file_len(raw_file)
	nstars=int( (nlines-3)/num_rows_per_object )

	# to hold mag/magerr info
	mag_arr=np.zeros(n_images)
	magerr_arr=np.zeros(n_images)


	# read in artstar data
	with open(raw_file,"r") as f:

		# Pass over header
		for x in range(3): f.readline()

		star_num=1
		while star_num <= nstars:
		
			# put all data for a star into a long list
			parts=list()
			for i in range(num_rows_per_object):
			        parts.extend(f.readline().split())
			
			
			# loop through star values and record them
			part_num=0
			
			id_el=int(parts[part_num])
			part_num+=1
			x_el=float(parts[part_num])
			part_num+=1
			y_el=float(parts[part_num])
			part_num+=1

			p_mag_magerr_pairs=int( (len(parts)-3-2)/2. )  # exclude id,x,y and chi,sharp
			for col in range(0, p_mag_magerr_pairs):
				mag_arr[col]=float(parts[part_num])
				part_num+=1
				magerr_arr[col]=float(parts[part_num])
				part_num+=1

			# now do chi, sharp
			chi_el=float(parts[part_num])
			part_num+=1
			sharp_el=float(parts[part_num])
			part_num+=1

			

			# see if star is an artstar by matching by coord
			start_loc=bisect.bisect_right(mag1_x_artstar,x_el-2) # start matching a little ahead of possible match


			found_match_w_artstar=0
			match_dist=1.0 # pix
			for el in range(start_loc,len(mag1_x_artstar)):
				# calc dist
				dx=mag1_x_artstar[el]-x_el
				dy=mag1_y_artstar[el]-y_el
				distance=np.sqrt(dx**2.+dy**2.)
				if distance <= match_dist:
					found_match_w_artstar=1
	
					artstar_mag1_input = mag1_artstar[el]
					artstar_mag2_input = mag2_artstar[el]

					break
				# end search if well-beyond possible match radius
				if mag1_x_artstar[el]-x_el > 2.0: break
			if found_match_w_artstar == 0: 
				star_num+=1
				continue	

			# correct for exptime differences
			for mag_el in range(n_images):
				mag_arr[mag_el] += exptime_mag_diff[mag_el]


			# convert to fluxes, get average mag
			avg_mag_band1,avg_magerr_band1 = get_avg_mag_magerr(mag_arr[0:second_band_index],magerr_arr[0:second_band_index])
			avg_mag_band2,avg_magerr_band2 = get_avg_mag_magerr(mag_arr[second_band_index:n_images],magerr_arr[second_band_index:n_images])

			star_num+=1




			# write matched info to file
			with open('condensed_artstar.dat',"a") as outf:
        		                	
				# write to file
				id_string=str(id_el)
				x_string =("{0:.3f}").format(x_el)
				y_string =("{0:.3f}").format(y_el)
				mag1_input_string=("{0:.3f}").format(artstar_mag1_input)
				mag2_input_string=("{0:.3f}").format(artstar_mag2_input)
				mag1_string=("{0:.3f}").format(avg_mag_band1)
				magerr1_string=("{0:.3f}").format(avg_magerr_band1)
				mag2_string=("{0:.3f}").format(avg_mag_band2)
				magerr2_string=("{0:.3f}").format(avg_magerr_band2)
				chi_string=("{0:.2f}").format(chi_el)
				sharp_string=("{0:.3f}").format(sharp_el)
				
				outf.write("   "+id_string+" "+\
					x_string+" "+y_string+" "+\
					mag1_input_string+" "+mag2_input_string+" "+\
					mag1_string+" "+magerr1_string+"   "+\
					mag2_string+" "+magerr2_string+"   "+\
					chi_string+"    "+sharp_string+"\n")




##################################################################


# User input
n_f814w_images=4
n_images=file_len("name.list")
second_band_index=n_f814w_images



# Write header to output file
os.system('rm -f condensed_artstar.dat')
with open('condensed_artstar.dat',"w") as f:
	f.write('# id,x,y,f814w_input,f606w_input,f814w_out,f814werr_out,f606w_out,f606werr_out,chi_out,sharp_out\n')



# Read filenames into list
filebasenames=list()
with open('name.list') as f:
	for row in f:
		# 'split' the string to aviod the newline \n form
		name=row.split()
		filebasenames.append(name[0])	


# Relevant files from which to extract data
f814w_add_files=glob.glob('*/*'+filebasenames[0]+'*.add')
f606w_add_files=glob.glob('*/*'+filebasenames[n_f814w_images-n_images]+'*.add')
raw_files=glob.glob('*/*_alf.raw')

# get exposure time from FITS image        
hdulist = pyfits.open(filebasenames[0]+'.fits')
ref_f814w_exptime=hdulist[0].header['EXPTIME']
hdulist = pyfits.open(filebasenames[n_f814w_images-n_images]+'.fits')
ref_f606w_exptime=hdulist[0].header['EXPTIME']

# get exposure time correction list
exptime_mag_diff=np.zeros(n_images)
for count,filebasename in enumerate(filebasenames):
	hdulist = pyfits.open(filebasename+'.fits')
	curr_exptime=hdulist[0].header['EXPTIME']
	if count < second_band_index:
		exptime_mag_diff[count] = np.log10(curr_exptime/ref_f814w_exptime)
	else:
		exptime_mag_diff[count] = np.log10(curr_exptime/ref_f606w_exptime)	

	

if len(f814w_add_files) != len(raw_files) or len(f606w_add_files) != len(raw_files): 
	print('number of files do not match') 
	exit(1)



for el in range(len(f814w_add_files)):

	print(raw_files[el],f814w_add_files[el],el)

	load_pull_write_artstars(raw_files[el],f814w_add_files[el],f606w_add_files[el],filebasenames,exptime_mag_diff,el,n_images,ref_f814w_exptime,ref_f606w_exptime,n_f814w_images)


