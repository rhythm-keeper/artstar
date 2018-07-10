pro condense_artstar_output


spawn,'rm -f condensed_artstar.dat'
openw,lunout,'condensed_artstar.dat',/get_lun
printf,lunout,'# id x y f606w_input f814w_input f606w_out f606werr_out f814w_out f814werr_out chi_out sharp_out'
free_lun,lunout



;spawn,'rm -f condensed_artstar.dat'

n_artstar_max=2e3

readcol,'name.list',filebasenames,format='(A)'

f814w_add_files=file_search('*/*'+filebasenames[0]+'*.add')
f606w_add_files=file_search('*/*'+filebasenames[-2]+'*.add')
raw_files=file_search('*/*_alf.raw')
mch_file='*_alf.mch'
n_f814w_images=6

; get reference exposure time so the right S/N is given to each frame
hdr=headfits(filebasenames[0]+'.fits')
ref_f814w_exptime=sxpar(hdr,'EXPTIME')
hdr=headfits(filebasenames[-2]+'.fits')
ref_f606w_exptime=sxpar(hdr,'EXPTIME')


if n_elements(f814w_add_files) ne n_elements(raw_files) or n_elements(f606w_add_files) ne n_elements(raw_files) then begin
	print,'number of files do not match' & stop
endif





for iter=0,n_elements(f814w_add_files)-1 do begin


	
	print,raw_files[iter],f814w_add_files[iter],iter

	load_pull_write_artstars,raw_files[iter],f814w_add_files[iter],f606w_add_files[iter],filebasenames,iter,mch_file,n_artstar_max,ref_f814w_exptime,ref_f606w_exptime,n_f814w_images


	

endfor ; iteration loop




end
