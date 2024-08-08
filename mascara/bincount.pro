function bincount, in_arr, weights=weights
;;; A substitute for numpy's bincount() function.
;;; Created 2024 Jul 23 by Daniel J. Stevens (UM-Duluth)
	upper = max(in_arr)
	counts = uintarr(upper + 1)
	print, n_elements(counts), upper
	for i=0, upper do begin
		if ~keyword_set(weights) then counts(i) = n_elements(where(in_arr eq i, /null)) $
		else begin
			mtch = where(in_arr eq i, /null)
			if mtch ne !NULL then counts[i] = total(weights[mtch]) else counts[i] = 0;counts[i] = total(weights[where(in_arr eq i)])
		endelse
	endfor
	return, counts
end
