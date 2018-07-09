pro gen_addstar_lf


; Image dimensions
max_x=4096
max_y=2048




TRGB_val=18.3 ; guess real TRGB val
F814W_cal=6.5 ; calibration 
M_TRGB=TRGB_val-F814W_cal ; TRGB val in instr. system for daophot

; resolution of magnitudes to input
nsamps=1./0.0005

rgb_slope=0.3
agb_slope=0.1


; RGB

; evenly spaced full mag range later to be sampled from
mags=findgen(nsamps)/(nsamps-1)*(M_TRGB+1-M_TRGB)+M_TRGB

rgb_lum_fcn=gen_lf(mags,nsamps,rgb_slope)

; AGB

mags=findgen(nsamps)/(nsamps-1)*(M_TRGB+1-(M_TRGB-1))+(M_TRGB-1)

agb_lum_fcn=gen_lf(mags,nsamps,agb_slope)

lum_fcn=rgb_lum_fcn


; Visualize the distribution of mags via histogram
;lum_hist,rgb_lum_fcn
;lum_hist,agb_lum_fcn,/oplot


; Star with 1500 RGB stars
rgb_subset=randomu(seed,1500)*n_elements(rgb_lum_fcn)
rgb_list=rgb_lum_fcn[rgb_subset]


; Want ratio of 4:1 at the tip
; Get current estimate of RGB
; Then add AGB luminosity functions
rgb_ind=where( abs(M_TRGB-rgb_list) lt 0.25 , num_rgb_at_TRGB)
current_ratio_at_tip=99.
agb_list=[]
agb_el=0L
star_incr=5
while current_ratio_at_tip gt 4 do begin

	; add agb stars in increments
	agb_list=[agb_list,agb_lum_fcn[agb_el:agb_el+star_incr-1]]

	agb_el+=star_incr	
	
	agb_ind=where( abs(M_TRGB-agb_list) lt 0.25 , num_agb_at_TRGB)

	current_ratio_at_tip=float(num_rgb_at_TRGB)/num_agb_at_TRGB
	;print,current_ratio_at_tip

endwhile

;lum_hist,rgb_list
;lum_hist,agb_list,/oplot

ninput=n_elements(rgb_list)+n_elements(agb_list)
n_rgb=n_elements(rgb_list)
n_agb=n_elements(agb_list)


; extra pixel range to span total width of image dithering
x_pix=randomu(seed,ninput)*(max_x+0)-0 ; 
y_pix=randomu(seed,ninput)*(max_y+0)-0 ; 
;x_pix=randomu(seed,ninput)*(max_x)
;y_pix=randomu(seed,ninput)*(max_y)



n_f814w_files=6

readcol,'name.list',filenamebase,format='(A)'

for i=0,n_elements(filenamebase)-1 do begin
	str="head -3 "+filenamebase[i]+".als  > "+filenamebase[i]+'.add'
	spawn,str
endfor


; read in distortion coefficients
; THIS FILE HAS BEEN HARD COPIED
; SINCE THE ALS FILENAMES HAVE THE
; ADDSTAR "A" APPPENDED TO THEIR NAMES
readcol,'1.mch',name_junk,junk2,AA_coeff,BB_coeff,CC_coeff,DD_coeff,EE_coeff,FF_coeff,$
	avg_offset,avg_offset_dispersion,$
	GG_coeff,HH_coeff,II_coeff,JJ_coeff,KK_coeff,LL_coeff,MM_coeff,NN_coeff,$
	OO_coeff,PP_coeff,QQ_coeff,RR_coeff,SS_coeff,TT_coeff,$
	format='(A,A,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D)'


; get reference exposure time so the right S/N is given to each frame
hdr=headfits(filenamebase[0]+'.fits')
ref_f814w_exptime=sxpar(hdr,'EXPTIME')
hdr=headfits(filenamebase[-2]+'.fits')
ref_f606w_exptime=sxpar(hdr,'EXPTIME')


; define the star colors
center_color=1.25
color_range=0.5
star_colors=center_color+(randomu(seed,n_rgb+n_agb)-0.5)*color_range


for file_loop=0,n_elements(filenamebase)-1 do begin

	
	; Adjust each star mag for exptime differences
		
	hdr=headfits(filenamebase[file_loop]+'.fits')
	curr_exptime=sxpar(hdr,'EXPTIME')
	if file_loop lt n_f814w_files then begin
		exptime_mag_diff=alog10(curr_exptime/ref_f814w_exptime) ; longer exp, brighter stars, subtract off this value
	endif else begin
		exptime_mag_diff=alog10(curr_exptime/ref_f606w_exptime)
	endelse


	openw,lun,filenamebase[file_loop]+'.add',/get_lun,/append

	for i=0,(n_rgb+n_agb)-1 do begin	

		; 6 order
		;common broyfunc_input,X0,Y0,A,B,C,D,E,F
		; 12 order
		;common broyfunc_input,X0,Y0,A,B,C,D,E,F,GG,HH,II,JJ,KK,LL
		; 20 order
		common broyfunc_input,X0,Y0,AA,BB,CC,DD,EE,FF,GG,HH,II,JJ,KK,LL,MM,NN,OO,PP,QQ,RR,SS,TT
		X0=x_pix[i]
		Y0=y_pix[i]
		AA=AA_coeff[file_loop]
		BB=BB_coeff[file_loop]
		CC=CC_coeff[file_loop]
		DD=DD_coeff[file_loop]
		EE=EE_coeff[file_loop]
		FF=FF_coeff[file_loop]
		GG=GG_coeff[file_loop]
		HH=HH_coeff[file_loop]
		II=II_coeff[file_loop]
		JJ=JJ_coeff[file_loop]
		KK=KK_coeff[file_loop]
		LL=LL_coeff[file_loop]
		MM=MM_coeff[file_loop]
		NN=NN_coeff[file_loop]
		OO=OO_coeff[file_loop]
		PP=PP_coeff[file_loop]
		QQ=QQ_coeff[file_loop]
		RR=RR_coeff[file_loop]
		SS=SS_coeff[file_loop]
		TT=TT_coeff[file_loop]

		; initial guess; just use the X,Y off values
		COORD = [x_pix[i]-AA,y_pix[i]-BB]
		
		; don't attempt transformation for the reference frame
		if file_loop ne 0 then begin
			new_coord = BROYDEN(COORD, 'BROYFUNC')
			x_pix_new=new_coord[0]
			y_pix_new=new_coord[1]
		endif else begin
			x_pix_new=x_pix[i]
			y_pix_new=y_pix[i]
		endelse


		if i lt n_rgb then begin
			; RGB stars; if lt n_f814w_files then use F814W else use F606W 1 mag fainter
			if file_loop lt n_f814w_files then curr_mag=rgb_list[i] else curr_mag=rgb_list[i]+star_colors[i]
                endif else begin
                        ; AGB stars 
                        if file_loop lt n_f814w_files then curr_mag=agb_list[i-n_rgb] else curr_mag=agb_list[i-n_rgb]+star_colors[i-n_rgb]
                endelse


		; add exptime adjustment
		curr_mag -= exptime_mag_diff


		if i lt n_rgb then id=long(i+2e5) else id=long(i+3e5-n_rgb)


		printf,lun,"   "+strtrim(id,2)+$
        	        " "+strtrim(string(x_pix_new,format='(F15.3)'),2)+$
        	        " "+strtrim(string(y_pix_new,format='(F15.3)'),2)+$
        	        "   "+strtrim(string(curr_mag,format='(F15.3)'),2)
	
	endfor

	free_lun,lun


endfor

end


FUNCTION broyfunc, COORD  
	common broyfunc_input
	
	
	;XY=( 2.*(COORD[0]-1.)/(RCOL-1.)-1. )*( 2.*(COORD[1]-1.)/(RROW-1.)-1. )
	;X2=1.5*( 2.*(COORD[0]-1.)/(RCOL-1.)-1.  )^2.-0.5
	;Y2=1.5*YS**2-0.5

	; Original linear transform
	;RETURN, [ -X0 + A + C*COORD[0] + E*COORD[1],$
        ;          -Y0 + B + D*COORD[0] + F*COORD[1] ]


	; 12 order transform, simplified
	;RCOL=4096
	;RROW=2048
        ;XS=2.*(COORD[0]-1.)/(RCOL-1.)-1.
        ;YS=2.*(COORD[1]-1.)/(RROW-1.)-1.
        ;XY=XS*YS
        ;X2=1.5*XS^2.-0.5
        ;Y2=1.5*YS^2.-0.5
	;RETURN, [ -X0 + A + C*COORD[0] + E*COORD[1] + GG*X2 + II*XY + KK*Y2,$
	;          -Y0 + B + D*COORD[0] + F*COORD[1] + HH*X2 + JJ*XY + LL*Y2 ]


	; 20 order transform, simplified

	RCOL=4096
        RROW=2048
        XS=2.*(COORD[0]-1.)/(RCOL-1.)-1.
        YS=2.*(COORD[1]-1.)/(RROW-1.)-1.
        XY=XS*YS
        X2=1.5*XS^2.-0.5
        Y2=1.5*YS^2.-0.5

	; 12 order transform, long form
	;RETURN, [ -X0 + A + C*COORD[0] + E*COORD[1] + GG*( 1.5*(2.*(COORD[0]-1.)/(RCOL-1.)-1.)^2.-0.5 ) + II*((2.*(COORD[0]-1.)/(RCOL-1.)-1.)*(2.*(COORD[1]-1.)/(RROW-1.)-1.)) + KK*(1.5*(2.*(COORD[1]-1.)/(RROW-1.)-1.)^2.-0.5),$	
        ;          -Y0 + B + D*COORD[0] + F*COORD[1] + HH*( 1.5*(2.*(COORD[0]-1.)/(RCOL-1.)-1.)^2.-0.5 ) + JJ*((2.*(COORD[0]-1.)/(RCOL-1.)-1.)*(2.*(COORD[1]-1.)/(RROW-1.)-1.)) + LL*(1.5*(2.*(COORD[1]-1.)/(RROW-1.)-1.)^2.-0.5) ]


	
	LINEAR_TERMS_X = AA + CC*COORD[0] + EE*COORD[1]
	LINEAR_TERMS_Y = BB + DD*COORD[0] + FF*COORD[1]
	QUAD_TERMS_X = GG*X2 + II*XY + KK*Y2
	QUAD_TERMS_Y = HH*X2 + JJ*XY + LL*Y2
	CUBIC_TERMS_X = ( MM*XS + OO*YS )*X2 + ( QQ*XS + SS*YS )*Y2
	CUBIC_TERMS_Y = ( NN*XS + PP*YS )*X2 + ( RR*XS + TT*YS )*Y2 

	RETURN, [-X0 + LINEAR_TERMS_X + QUAD_TERMS_X + CUBIC_TERMS_X,$
		 -Y0 + LINEAR_TERMS_Y + QUAD_TERMS_Y + CUBIC_TERMS_Y ]
        ;RETURN, [ -X0 + A + C*COORD[0] + E*COORD[1] + GG*X2 + II*XY + KK*Y2 + ( MM*XS + OO*YS )*X2 + ( QQ*XS + SS*YS )*Y2 ,$
        ;          -Y0 + B + D*COORD[0] + F*COORD[1] + HH*X2 + JJ*XY + LL*Y2 + ( NN*XS + PP*YS )*X2 + ( RR*XS + TT*YS )*Y2 ]



END  
