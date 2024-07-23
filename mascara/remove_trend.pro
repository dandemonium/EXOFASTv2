function remove_trend, time, flux, flux_err, lstseq, lst, x, y, sky, trend, mask, $
												nobs, modelflux
;;; Port of the MASCARA team's remove_trend.py to IDL.
;;;
;;; INPUTS:
;;;   time, flux, flux_err: usual light curve params
;;;   lstseq, lst, x, y, sky, loclin: MASCARA detrending parameters (floats)
;;;   mask: boolean MASCARA detrending parameter.
;;;   nobs: number of points in the light curve.
;;;   modelflux: transit model flux.
	strides = [0, cumsum(nobs)]
    
    for i = 0, n_elements(nobs)-1 do begin
        i1 = strides[i]
        i2 = strides[i+1]
        ; check for data in this block. Note: remove_trend.py has "if np.all(~mask): continue"
		; which is only True if *all* mask elements are False, but I think we just want *some*
		; unmasked data in this range?
        mask_new = mask[i1:i2]       
        if array_equal(mask_new, 0) eq 1 then continue ; following remove_trend.py verbatim for now.
        
        tmptime = time[i1:i2][where(mask_new)]
        tmpflux = flux[i1:i2][where(mask_new)]
        tmperr = flux_err[i1:i2][where(mask_new)]
        
        if n_elements(model) eq 0 then flux = tmpflux else flux = tmpflux - modelflux[i1:i2][where(mask_new)]
        trend, newmask, chi2 = mascara_loclin, time, flux, flux_err, lstseq, x, y, sky, window=5., maxiter=50, dtol=1e-3        
        trend[i1:i2][where(mask_new)] = trend
    endfor
	return, [newmask, chi2]
end