pro load_pull_write_artstars,raw_file,f814w_add_file,f606w_add_file,filebasenames,iter,mch_file,n_artstar_max,ref_f814w_exptime,ref_f606w_exptime,n_f814w_images


	; search coordinate list of input artificial stars
	readcol,f814w_add_file,id_artstar,x_artstar,y_artstar,skipline=3


	; make monotonically increasing
	x_sorted=sort(x_artstar)
	id_artstar=id_artstar[x_sorted]
	x_artstar=x_artstar[x_sorted]
	y_artstar=y_artstar[x_sorted]

;ind=where(id_artstar eq 300089)
;print,ind
;print,x_artstar[ind],y_artstar[ind]
;stop




        n_images=file_lines(mch_file)

        openr,lun,raw_file,/get_lun
        linenum=0L
        line=''
        header=[]
        while linenum lt 3 do begin
                linenum+=1
                readf,lun,line,format='(A)'
                header=[header,line]
        endwhile
        delvar,header

        f814w_out=fltarr(n_artstar_max)-99
        f814werr_out=fltarr(n_artstar_max)-99
        f606w_out=fltarr(n_artstar_max)-99
        f606werr_out=fltarr(n_artstar_max)-99
        id_out=fltarr(n_artstar_max)-99
        x_out=fltarr(n_artstar_max)-99
        y_out=fltarr(n_artstar_max)-99
        chi_out=fltarr(n_artstar_max)-99
        sharp_out=fltarr(n_artstar_max)-99
        mag_out_el=0L


        while ~eof(lun) do begin


                readf,lun,line,format='(A)'

                id_el=float(gettok(line," "))
                x_el=float(gettok(line," "))
                y_el=float(gettok(line," "))


                mag=[]
                magerr=[]



                curr_line=0
                while curr_line lt n_images do begin

                        ; six value is set by the number of 
                        if curr_line mod 6 eq 0 and curr_line ne 0 then readf,lun,line

                        mag=[mag,float(gettok(line," "))]
                        magerr=[magerr,float(gettok(line," "))]

                        curr_line+=1

                endwhile


                ;readf,lun,line

                chi_el=float(gettok(line," "))
                sharp_el=float(gettok(line," "))

		;print,id_el,x_el,y_el
		;print,mag
		;print,'------'
		;
		;if id_el eq 503 then stop

		
		
		start_loc=value_locate(x_artstar,x_el)
		start_loc=start_loc[0] ; for some reason needs this for loop to work below...
		; boundary conditions http://www.harrisgeospatial.com/docs/VALUE_LOCATE.html
		if start_loc eq -1 then start_loc=0


		found_match_w_artstar=0
		match_dist=1.0 ; pixel
		; loop thourgh artstars
		;print,'loop through artstars'
		for el=start_loc,n_elements(x_artstar)-1 do begin
			; calc dist
			dx=x_artstar[el]-x_el
			dy=y_artstar[el]-y_el
			distance=sqrt(dx^2.+dy^2.)
			min_dist=min(distance)
		;	if id_el eq 503 then begin
		;
		;		print,x_el,y_el,min_dist,id_el
		;		print,el
		;		print,x_artstar[180],y_artstar[180]
		;	endif
			if min_dist le match_dist then begin
				found_match_w_artstar=1
				break	
			endif
			if x_artstar[el]-x_el gt 2.0 then break ; stop searching beyond some fwhm of stars
		endfor
		;print,'end loop through artstars'
		if found_match_w_artstar eq 0 then continue

	
                ;if id_el lt 2e5 or id_el ge 5e5 then continue ; non valid IDs


                ; Average magnitudes F814W
                mag_list=[]
                magerr_list=[]
                for mag_loop=0,n_f814w_images-1 do begin

                  if mag[mag_loop] ge 99. or mag[mag_loop] lt 9 then continue

                  hdr=headfits(filebasenames[mag_loop]+'.fits')
                  curr_exptime=sxpar(hdr,'EXPTIME')
                  exptime_mag_diff=alog10(curr_exptime/ref_f814w_exptime) ; add back on this value

                  mag_list=[mag_list,mag[mag_loop] + exptime_mag_diff ]
                  magerr_list=[magerr_list,magerr[mag_loop] ]

                endfor ; mag loop


                if n_elements(mag_list) eq 0 then continue

                ; average by flux
                fluxes=[]
                fluxerrs=[]
                for mag_loop=0,n_elements(mag_list)-1 do begin
                        flux=10.^( (mag_list[mag_loop]-25.)/(-2.5) )
                        fluxes=[fluxes,flux]
                        fluxerrs=[fluxerrs,flux*magerr_list[mag_loop]/1.0857]
                endfor
                ; remove outliers 
                ;mmad=mad(fluxes)
                ;remove_outliers=where( abs(fluxes - median(fluxes,/even)) le 2.0*mmad )
                ;fluxes=fluxes[remove_outliers]
                ;fluxerrs=fluxerrs[remove_outliers]

                wflux=wmean(fluxes,fluxerrs)
                f814w_out[mag_out_el]=-2.5*alog10(wflux)+25
                wfluxerr=1./sqrt(total(1./fluxerrs^2.))
                wmagerr=wfluxerr*1.0857/wflux

                f814werr_out[mag_out_el]=wmagerr



                delvar,fluxes
                delvar,fluxerrs
                delvar,mag_list
                delvar,magerr_list



                ; Average magnitudes F606W
                mag_list=[]
                magerr_list=[]
                for mag_loop=n_f814w_images,n_images-1 do begin
                  if mag[mag_loop] ge 99. or mag[mag_loop] lt 9 then continue

                 hdr=headfits(filebasenames[mag_loop]+'.fits')
                 curr_exptime=sxpar(hdr,'EXPTIME')
                 exptime_mag_diff=alog10(curr_exptime/ref_f606w_exptime) ; add back on this value

                  mag_list=[mag_list,mag[mag_loop] + exptime_mag_diff ]
                  magerr_list=[magerr_list,magerr[mag_loop] ]
                endfor ; mag loop

                if n_elements(mag_list) eq 0 then begin
                        mag_list=fltarr(n_images-n_f814w_images)-99
                        magerr_list=fltarr(n_images-n_f814w_images)-99
                endif

                ; average by flux
                fluxes=[]
                fluxerrs=[]
                for mag_loop=0,n_elements(mag_list)-1 do begin
                        flux=10.^( (mag_list[mag_loop]-25.)/(-2.5) )
                        fluxes=[fluxes,flux]
                        fluxerrs=[fluxerrs,flux*magerr_list[mag_loop]/1.0857]
                endfor
                ; remove outliers 
                ;mmad=mad(fluxes)
                ;remove_outliers=where( abs(fluxes - median(fluxes,/even)) le 2.0*mmad )
                ;fluxes=fluxes[remove_outliers]
                ;fluxerrs=fluxerrs[remove_outliers]

                ;wflux=wmean(fluxes,fluxerrs)
                ;f606w_out[mag_out_el]=-2.5*alog10(wflux)+25
                ;wmagerr=stdev(fluxerrs)/wflux*1.0857/sqrt(n_elements(fluxes)-1)        
                wflux=wmean(fluxes,fluxerrs)
                f606w_out[mag_out_el]=-2.5*alog10(wflux)+25
                wfluxerr=1./sqrt(total(1./fluxerrs^2.))
                wmagerr=wfluxerr*1.0857/wflux


                f606werr_out[mag_out_el]=wmagerr


                chi_out[mag_out_el]=chi_el
                sharp_out[mag_out_el]=sharp_el
                id_out[mag_out_el]=id_el
                x_out[mag_out_el]=x_el
                y_out[mag_out_el]=y_el
                mag_out_el+=1


                delvar,fluxes
                delvar,fluxerrs
                delvar,mag_list
                delvar,magerr_list


                delvar,mag
                delvar,magerr

        endwhile

	close,lun
        free_lun,lun

        ind=where(f814w_out ne -99)
        f814w_out=f814w_out[ind]
        f814werr_out=f814werr_out[ind]
        f606w_out=f606w_out[ind]
        f606werr_out=f606werr_out[ind]
        chi_out=chi_out[ind]
        sharp_out=sharp_out[ind]
        id_out=id_out[ind]
        x_out=x_out[ind]
        y_out=y_out[ind]


	; sort by x-coord for fast searching below
	ind=sort(x_out)
	f814w_out=f814w_out[ind]
        f814werr_out=f814werr_out[ind]
        f606w_out=f606w_out[ind]
        f606werr_out=f606werr_out[ind]
        chi_out=chi_out[ind]
        sharp_out=sharp_out[ind]
        id_out=id_out[ind]
        x_out=x_out[ind]
        y_out=y_out[ind]


        readcol,f606w_add_file,id_in,x_in,y_in,f606w_in,skipline=3
        readcol,f814w_add_file,id_in,x_in,y_in,f814w_in,skipline=3



	print,'matching raw mags with artstars input'

	openw,lunout,'condensed_artstar.dat',/get_lun,/append

        for i=0,n_elements(f814w_in)-1 do begin


		;dx=x_out-x_in[i]
		;dy=y_out-y_in[i]
		;distances=sqrt(dx^2.+dy^2.)
		;min_dist=min(distances)
		;ind=where(distances eq min_dist)
		;if min_dist le 1 then matched=1 else matched=0


                ;ind=where(id_in[i] eq id_out,matched)

                ;all_id=id_in[i]
                ;all_f814w_in=f814w_in[i]
                ;all_f606w_in=f606w_in[i]

                ;if id_in eq 200000 then begin
                ;        print,f814w_in[i],f814w_out[ind]
                ;       stop
                ;endif

		start_loc=value_locate(x_out,x_in[i])
                start_loc=start_loc[0] ; for some reason needs this for loop to work below...
                ; boundary conditions http://www.harrisgeospatial.com/docs/VALUE_LOCATE.html
                if start_loc eq -1 then start_loc=0

		;print,start_loc,n_elements(x_out)
                found_match_w_artstar=0
                ; loop thourgh artstars
                for el=start_loc,n_elements(x_out)-1 do begin
                        ; calc dist
                        dx=x_out[el]-x_in[i]  
                        dy=y_out[el]-y_in[i]
                        distance=sqrt(dx^2.+dy^2.)
                        min_dist=min(distance)
                        if min_dist le match_dist then begin
                                found_match_w_artstar=1
				ind=el
                                break
                        endif
                        if x_out[el]-x_el gt 2.0 then break ; stop searching beyond some fwhm of stars
                endfor



                if found_match_w_artstar then begin
                        all_f814w_out=f814w_out[ind]
                        all_f814werr_out=f814werr_out[ind]
                        all_f606w_out=f606w_out[ind]
                        all_f606werr_out=f606werr_out[ind]
                        all_chi_out=chi_out[ind]
                        all_sharp_out=sharp_out[ind]
                        ;all_x_out=x_out[ind]
                        ;all_y_out=y_out[ind]
                endif else begin
                        all_f814w_out=-999
                        all_f814werr_out=-999
                        all_f606w_out=-999
                        all_f606werr_out=-999
                        all_chi_out=-999
                        all_sharp_out=-999


                endelse




                printf,lunout,strtrim(id_in[i],2)+" "+$
                        strtrim(x_in[i],2)+" "+$
                        strtrim(y_in[i],2)+" "+$
                        strtrim(f606w_in[i],2)+" "+$
                        strtrim(f814w_in[i],2)+" "+$
                        strtrim(all_f606w_out,2)+" "+$
                        strtrim(all_f606werr_out,2)+" "+$
                        strtrim(all_f814w_out,2)+" "+$
                        strtrim(all_f814werr_out,2)+" "+$
                        strtrim(all_chi_out,2)+" "+$
                        strtrim(all_sharp_out,2)


        endfor

	print,'end matching raw mags with artstars input'

	close,lunout
        free_lun,lunout


        delvar,f814w_out
        delvar,f814werr_out
        delvar,f606w_out
        delvar,f606werr_out
        delvar,id_out
        delvar,x_out
        delvar,y_out
        delvar,chi_out
        delvar,sharp_out





;return


end
