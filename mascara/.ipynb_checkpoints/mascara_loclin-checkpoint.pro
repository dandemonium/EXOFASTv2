;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Detrending functions.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
FUNCTION LOCAL_LINEAR, jd, lstseq, x, y, sky, mag, emag, window=window, maxiter=maxiter, dtol=dtol
    if ~keyword_set(window) then window=5d0
	if ~keyword_set(maxiter) then maxiter=50
	if ~keyword_set(dtol) then dtol=1e-3
    lstidx = (lstseq MOD 270)
    
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
	end
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
	end
FUNCTION linfit, lstidx, x, y, sky, mag, emag
    sorted = sort(lstidx)
    invsort = sort(sorted)
	
    lstidx = lstidx[sorted]
    x = x[sorted]
    y = y[sorted]
    sky = sky[sorted]
    mag = mag[sorted]
    emag = emag[sorted]
    
    idx = unique(lstidx)
    
    nobs = pdf(idx)
    strides = append(0, total(nobs, /cumulative))
    
    xbar = BincountWithWeights(idx, x)/Bincount(idx)
    ybar = BincountWithWeights(idx, y)/Bincount(idx)
    
    mat = np.column_stack([np.ones_like(mag), x-xbar[idx], y-ybar[idx], sky])
    
    pars = fltarr((n_elements(nobs), 4))
    pars[:,0] = BincountWithWeights(idx, mag/emag**2)/BincountWithWeights(idx, 1/emag**2)
    
    for i=0, n_elements(len(nobs)) do begin        
        if nobs[i] lt 5 then continue            
        i1 = strides[i]
        i2 = strides[i+1]
        
        pars[i] = np.linalg.lstsq(mat[i1:i2]/emag[i1:i2,np.newaxis], mag[i1:i2]/emag[i1:i2], rcond=None)[0]
	endfor
    trend = np.sum(pars[idx]*mat, axis=1)

    return, trend[invsort], (nobs > 4)[idx][invsort] 
	end
FUNCTION Bincount, inArray, inMinLength=1
    COMPILE_OPT idl2

    ; Ensure inArray is an integer type
    IF (SIZE(inArray, /TNAME) NE 'LONG' AND SIZE(inArray, /TNAME) NE 'ULONG' AND $
        SIZE(inArray, /TNAME) NE 'INT' AND SIZE(inArray, /TNAME) NE 'UINT' AND $
        SIZE(inArray, /TNAME) NE 'BYTE') THEN $
        MESSAGE, 'Input array must contain integers', /CONTINUE

    maxValue = MAX(inArray)
    IF (maxValue LT 0) THEN BEGIN
        ; No positive values, return an empty array
        outArray = LONARR(0)
        RETURN
    ENDIF

    outArraySize = MAX([maxValue + 1, inMinLength])
    clippedArray = CLIP(inArray, 0, maxValue)

    outArray = LONARR(outArraySize)
    FOR i = 0, N_ELEMENTS(clippedArray) - 1 DO outArray[clippedArray[i]]++
	return, outArray
	end
FUNCTION BincountWithWeights, inArray, inWeights, inMinLength=1
    COMPILE_OPT idl2

    ; Ensure inArray is an integer type
    IF (SIZE(inArray, /TNAME) NE 'LONG' AND SIZE(inArray, /TNAME) NE 'ULONG' AND $
        SIZE(inArray, /TNAME) NE 'INT' AND SIZE(inArray, /TNAME) NE 'UINT' AND $
        SIZE(inArray, /TNAME) NE 'BYTE') THEN $
        MESSAGE, 'Input array must contain integers', /CONTINUE

    ; Ensure inArray and inWeights have the same shape
    IF (SIZE(inArray, /DIMENSIONS) NE SIZE(inWeights, /DIMENSIONS)) THEN $
        MESSAGE, 'Weights array must be the same shape as the input array.', /CONTINUE

    maxValue = MAX(inArray)
    IF (maxValue LT 0) THEN BEGIN
        ; No positive values, return an empty array
        outArray = LONARR(0)
        RETURN
    ENDIF

    outArraySize = MAX([maxValue + 1, inMinLength])
    clippedArray = CLIP(inArray, 0, maxValue)

    outArray = LONARR(outArraySize)
    FOR i = 0, N_ELEMENTS(clippedArray) - 1 DO outArray[clippedArray[i]] += inWeights[i]
	return, outArray
	end

;;; move below to exofastv2?
function REMOVE_TREND, lc, nobs, model, method='loclin', options
    strides = [0, CUMSUM(nobs)]
    
    FOR i = 0, N_ELEMENTS(nobs)-1 DO BEGIN
        i1 = strides[i]
        i2 = strides[i+1]
        
        mask = lc.mask[i1:i2]
        
        IF TOTAL(mask) EQ 0 THEN CONTINUE
        
        tmp = lc[i1:i2][WHERE(mask)]
        
        IF N_ELEMENTS(model) EQ 0 THEN mag = tmp.mag ELSE mag = tmp.mag - model[i1:i2][WHERE(mask)]
        trend, flag, chisq = LOCAL_LINEAR, jd, lstseq, x, y, sky, mag, emag, window=5., maxiter=50, dtol=1e-3
        
        lc.trend[i1:i2][WHERE(mask)] = trend
    ENDFOR
END
