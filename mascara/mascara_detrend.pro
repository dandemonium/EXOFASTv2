PRO WRAP_LST, lst, result
    ; Wrap the LST
    sort = SORT(lst)
    gap = MAX(DIFF(lst[sort]))
    arg = WHERE(DIFF(lst[sort]) EQ gap, cnt)
    gap0 = 24.0 - (MAX(lst) - MIN(lst)) ; Gap across lst=0.

    IF (gap GT gap0) THEN BEGIN
        lst = (lst - lst[sort[arg[0]+1]]) MOD 24.0
    ENDIF
    
    result = lst
END

PRO SCALE, x, result
    ; Scale a coordinate to run from -1 to 1
    xmax = MAX(x)
    xmin = MIN(x)
    x = (x - xmin) / (xmax - xmin)
    x = 2.0 * x - 1.0
    result = x
END

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Detrending methods.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

FUNCTION POLYNOMIAL, jd, mag, emag, deg, trend, mat, pars, chisq
    ; Default degree
    IF N_ELEMENTS(deg) EQ 0 THEN deg = 0

    ; Create the matrix
    SCALE, jd, x
    mat = LEGAUSS(x, deg)
    
    ; Compute the best fit
    A = mat / REFORM(emag, 1, N_ELEMENTS(emag))
    B = mag / emag
    pars = LINFIT(A, B)
    
    ; Evaluate the best fit
    trend = TOTAL(pars * mat, 2)
    chisq = TOTAL(((mag - trend) / emag)^2)
END

FUNCTION LEGENDRE_DETREND, jd, lst, sky, mag, emag, s_jd, s_lst, sig, maxiter, trend, mat, pars, chisq
    ; Default parameters
    IF N_ELEMENTS(s_jd) EQ 0 THEN s_jd = 5.0
    IF N_ELEMENTS(s_lst) EQ 0 THEN s_lst = 0.25
    IF N_ELEMENTS(sig) EQ 0 THEN sig = 10.0
    IF N_ELEMENTS(maxiter) EQ 0 THEN maxiter = 50

    ; Long-term variations
    deg0 = FLOOR((MAX(jd) - MIN(jd)) / s_jd)
    SCALE, jd, x
    mat0 = LEGAUSS(x, deg0)[*,1:*]
    
    ; PSF variations
    WRAP_LST, lst, lst
    deg1 = FLOOR((MAX(lst) - MIN(lst)) / s_lst)
    SCALE, lst, x
    mat1 = LEGAUSS(x, deg1)
    
    ; Sky dependency
    sky = REFORM(sky, 1, N_ELEMENTS(sky))
    mat1_sky = mat1 * sky / MAX(sky)
    
    ; Combine matrices
    mat = [mat0, mat1, mat1_sky]
    
    ; Solve iteratively
    mask = BYTARR(N_ELEMENTS(jd), /BOOL)
    mask[0] = 1B
    
    FOR niter = 0, maxiter-1 DO BEGIN
        ; Compute the best fit
        A = mat[mask, *] / REFORM(emag[mask], 1, N_ELEMENTS(emag[mask]))
        B = mag[mask] / emag[mask]
        pars = LINFIT(A, B)
        fit = TOTAL(pars * mat, 2)
        
        ; Compute the residuals
        res = mag - fit
        m0 = MEDIAN(res)
        m1 = MEDABSDEV(res)
        
        ; Update the mask
        old_mask = mask
        mask = ABS(res - m0) LT sig * m1
        
        ; Check convergence
        IF TOTAL(mask EQ old_mask) EQ N_ELEMENTS(mask) THEN BREAK
        IF TOTAL(mask) EQ 0 THEN BREAK
    ENDFOR
    
    ; Evaluate the best fit
    trend = TOTAL(pars * mat, 2)
    chisq = TOTAL(((mag - trend) / emag)^2)
END

FUNCTION LOCAL_LINEAR, jd, lstseq, x, y, sky, mag, emag, window=5., maxiter=50, dtol=1e-3
    
    lstidx = (lstseq % 270)
    
    trend = replicate(0, n_elements(mag))
    trend0 = replicate(0, n_elements(mag))
    trend1 = replicate(0, n_elements(mag))
    
    for niter=0, maxiter do begin            
        trend0, mask = linfit(lstidx, x, y, sky, mag - trend1, emag)
        trend1 = moving_mean(jd, mag - trend0, emag, window)        
        if niter gt 0 then begin        
		    trend_diff = abs(trend - trend0 - trend1)
            if trend_diff < dtol then begin
			    print, trend_diff < dtol
                break
			endif
        endif    
        trend = trend0 + trend1
     endfor       
    ; Evaluate the best fit.  
    trend = trend0 + trend1
    chisq = total(((mag - trend)/emag)^2)
        
    return, trend, mask, chisq 

function moving_mean, x, y, yerr=yerr, window=3d0
    ; Compute a moving mean along the x-axis.

    ; Set the weights.
    if ~keyword_set(yerr) then weights = replicate(1, n_elements(y))
    else weights = 1/yerr**2

    ; Sums for computing the mean.
    sum1 = append(0, total(weights*y, /cumulative))
    sum2 = append(0, total(weights, /cumulative))

    ; Indices at the start and end of the window. Assumes time array
	; x is already sorted.
    i = where(x ge (x - window/2.))[0]
    j = where(x gt (x + window/2.))[0]

    ; Compute the mean.
    avg = (sum1[j] - sum1[i])/(sum2[j] - sum2[i])
    
    return, avg


PRO REMOVE_TREND, lc, nobs, model, method='loclin', options
    strides = [0, CUMSUM(nobs)]
    
    FOR i = 0, N_ELEMENTS(nobs)-1 DO BEGIN
        i1 = strides[i]
        i2 = strides[i+1]
        
        mask = lc.mask[i1:i2]
        
        IF TOTAL(mask) EQ 0 THEN CONTINUE
        
        tmp = lc[i1:i2][WHERE(mask)]
        
        IF N_ELEMENTS(model) EQ 0 THEN mag = tmp.mag ELSE mag = tmp.mag - model[i1:i2][WHERE(mask)]
        trend, flag, chisq = LOCAL_LINEAR, jd, lstseq, x, y, sky, mag, emag, window=5., maxiter=50, dtol=1e-3
;        CASE method OF
;            'none': 
;            'polynomial': POLYNOMIAL, tmp.jd, mag, tmp.emag, options.deg, trend, mat, pars, chisq
;            'legendre': LEGENDRE_DETREND, tmp.jd, tmp.lst, tmp.sky, mag, tmp.emag, options.s_jd, options.s_lst, options.sig, options.maxiter, trend, mat, pars, chisq
;            'loclin': LOCAL_LINEAR, jd, lstseq, x, y, sky, mag, emag, window=5., maxiter=50, dtol=1e-3
               
            ;'fourier': ; not supported
                ; Add FOURIER function implementation here
;            ELSE:
;                PRINT, 'Unknown detrending method "', method, '"'
;                RETURN
        ENDCASE
        
        lc.trend[i1:i2][WHERE(mask)] = trend
    ENDFOR
END
```
