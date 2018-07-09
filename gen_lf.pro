function gen_lf,mags,nsamps,slope


; the function that describes the
; abundance of stars vs mag is
; equivalent to a prob function
xvals=findgen(nsamps)/(nsamps-1)
prob_mass=10.^(slope*xvals)
; normalize to 1
prob_mass = prob_mass / total(prob_mass)

; To get mag distr that follows
; prob mass function, need to compute CDF 
cdf=fltarr(n_elements(prob_mass))
for el=0,n_elements(cdf)-1 do begin
        cdf[el]=total(prob_mass[0:el])
endfor
; Plot CDF
;cgplot,mags,cdf,psym=3
;stop

; randomly select mags from cdf
; via inverse transform sampling
nrands=200e3
rand_vals=randomu(seed,nrands)


lum_fcn=fltarr(nrands)
count=0L
for el=0,nrands-1 do begin

        diff_dist=abs(rand_vals[el]-cdf)
        ind=where(diff_dist eq min(diff_dist))

        if n_elements(ind) gt 1 then begin
                 shuffled_ind=shuffle(ind)
                 ind=ind[0]
        endif

        lum_fcn[count]=mags[ind]
        count+=1

endfor



return,lum_fcn




end
