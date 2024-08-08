function mascara_linfit, lstidx, x, y, sky, mag, emag
    sorted = sort(lstidx)
    invsort = sort(sorted)
	
    lstidx = lstidx[sorted]
    x = x[sorted]
    y = y[sorted]
    sky = sky[sorted]
    mag = mag[sorted]
    emag = emag[sorted]
    
	idx_inv = invert_uniq(lstidx)
	nobs = bincount(idx_inv)
    strides = append(0, total(nobs, /cumulative))
    
    xbar = bincount(idx_inv, weights=x)/bincount(idx_inv)
    ybar = bincount(idx_inv, weights=y)/bincount(idx_inv)
    mat = transpose([dblarr(n_elements(mag))], [x-xbar[idx_inv]], [y-ybar[idx_inv]], [sky]])
    mat[0,*] = 1d0
    pars = dblarr((n_elements(nobs), 4))
    pars[:,0] = bincount(idx_inv, weights=mag/emag**2)/bincount(idx_inv, weights=1/emag**2)
    
    for i=0, n_elements(len(nobs)) do begin        
        if nobs[i] lt 5 then continue            
        i1 = strides[i]
        i2 = strides[i+1]
;;; EVERYTHING BELOW THIS LINE NEEDS TO BE CONVERTED        
        pars[i] = LA_LEAST_SQUARES(mat[i1:i2]/emag[i1:i2,np.newaxis], mag[i1:i2]/emag[i1:i2], rcond=0)[0]
	endfor
    trend = total(pars[idx_inv]*mat, axis=1)

    return, trend[invsort], (nobs > 4)[idx_inv][invsort] 
end