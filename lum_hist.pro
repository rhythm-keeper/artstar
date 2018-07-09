pro lum_hist,lum_fcn,oplot=oplot

; Visualize the distribution of mags via histogram
min_val=min(lum_fcn)
max_val=max(lum_fcn)

binsize=0.025
nbins=long((max_val-min_val)/binsize )+1
x_pos=[min_val]
y_pos=[0]
for el=0,nbins-1 do begin
        ind=where( lum_fcn gt min_val+el*binsize and $
                 lum_fcn lt min_val+(el+1)*binsize ,count)

        x_pos=[x_pos,min_val+el*binsize]
        x_pos=[x_pos,min_val+(el+1)*binsize]

        y_pos=[y_pos,count]
        y_pos=[y_pos,count]
endfor

; last value going back to zero
x_pos=[x_pos,min_val+el*binsize]
y_pos=[y_pos,0]

if keyword_set(oplot) then cgoplot,x_pos,y_pos,color=cgcolor('black'),thick=2.0 $
else cgplot,x_pos,y_pos,color=cgcolor('black'),thick=2.0



end
