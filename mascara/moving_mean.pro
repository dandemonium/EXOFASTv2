function moving_mean, x, y, yerr=yerr, wndw=wndw
    ; Compute a moving mean along the x-axis.
	if ~keyword_set(wndw) then wndw=3d0
    ; Set the weights.
    if ~keyword_set(yerr) then weights = replicate(1, n_elements(y))
    else weights = 1/yerr**2

    ; Sums for computing the mean.
    sum1 = append(0, total(weights*y, /cumulative))
    sum2 = append(0, total(weights, /cumulative))

    ; Indices at the start and end of the wndw. Assumes time array
	; x is already sorted.
    i = where(x ge (x - wndw/2.))[0]
    j = where(x gt (x + wndw/2.))[0]

    ; Compute the mean.
    avg = (sum1[j] - sum1[i])/(sum2[j] - sum2[i])
    
    return, avg
	end