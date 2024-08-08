function mascara_loclin, time, flux, flux_err, lstseq, x, y, sky, wndw=wndw, maxiter=maxiter, dtol=dtol
;;; IDL port of the MASCARA team's Python code to perform "local linear" detrending.
;;; Inputs:
;;;   time, flux, flux_err: standard light curve inputs
;;;   lstseq, x, y, sky: MASCARA-specific detrending parameters
;;;   wndw=wndw: size of window for moving_mean detrending
;;;   maxiter=maxiter: maximum number of line-fitting iterations to perform
;;;   dtol=dtol: tolerance for line fitting.
;;; Written by: Daniel J. Stevens (UM-Duluth), 2024 July 23
	if ~keyword_set(wndw) then wndw=5d0
	if ~keyword_set(maxiter) then maxiter=50
	if ~keyword_set(dtol) then dtol=1e-3
    lstidx = (lstseq mod 270)
    
    trend = replicate(0, n_elements(flux))
    trend0 = replicate(0, n_elements(flux))
    trend1 = replicate(0, n_elements(flux))
	mag = -2.5*alog(flux)
	emag = fluxerr * 2.5/(flux * alog(10))
	
    for niter=0, maxiter do begin            
        trend0, mask = linfit(lstidx, x, y, sky, mag - trend1, emag)
        trend1 = moving_mean(time, mag - trend0, emag, wndw)              
        if niter gt 0 then begin        
		    trend_diff = abs(trend - trend0 - trend1)
            if trend_diff < dtol then begin
;			    print, trend_diff < dtol
                break
			endif
        endif    
        trend = trend0 + trend1
     endfor       
    ; Evaluate the best fit.  
    trend = trend0 + trend1
    chisq = total(((mag - trend)/emag)^2)
        
    return, trend, mask, chisq 
end

