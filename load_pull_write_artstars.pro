pro load_pull_write_artstars,raw_files,f814w_add_files,f606w_add_files,filebasenames,iter,mch_file,n_artstar_max,ref_f814w_exptime,ref_f606w_exptime,n_f814w_images



        n_images=file_lines(mch_file)

        openr,lun,raw_files[iter],/get_lun
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


		; search coordinate list of input artificial stars
		start_loc=value_locate(x_artstars,x_el)
		; loop thourgh artstars
		for el=start_loc,n_elements(x_artstars)-1 do begin
			; calc dist
			if x_artstars[el]-x_el gt 2*1.6 then el=n_elements(x_artstars) ; stop searching beyond some fwhm of stars
		endfor



                if id_el lt 2e5 or id_el ge 5e5 then continue ; non valid IDs


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

                ;if id_el eq 200000 then begin
                ;        print,mag_list 
                ;       print,f814w_out[mag_out_el]
                ;       stop
                ;endif



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




        readcol,f606w_add_files[iter],id_in,x_in,y_in,f606w_in,skipline=3
        readcol,f814w_add_files[iter],id_in,x_in,y_in,f814w_in,skipline=3



        for i=0,n_elements(f814w_in)-1 do begin
                ind=where(id_in[i] eq id_out,matched)

                all_id=id_in[i]
                all_f814w_in=f814w_in[i]
                all_f606w_in=f606w_in[i]

                ;if id_in eq 200000 then begin
                ;        print,f814w_in[i],f814w_out[ind]
                ;       stop
                ;endif


                if matched then begin
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

                all_x_out=x_in[i]
                all_y_out=y_in[i]


                openw,lunout,'condensed_artstar.dat',/get_lun,/append
                printf,lunout,strtrim(all_id,2)+" "+$
                        strtrim(all_x_out,2)+" "+$
                        strtrim(all_y_out,2)+" "+$
                        strtrim(all_f606w_in,2)+" "+$
                        strtrim(all_f814w_in,2)+" "+$
                        strtrim(all_f606w_out,2)+" "+$
                        strtrim(all_f606werr_out,2)+" "+$
                        strtrim(all_f814w_out,2)+" "+$
                        strtrim(all_f814werr_out,2)+" "+$
                        strtrim(all_chi_out,2)+" "+$
                        strtrim(all_sharp_out,2)
                free_lun,lunout




                ;count+=1

        endfor


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
