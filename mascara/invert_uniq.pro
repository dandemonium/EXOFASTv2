function invert_uniq, arr
;;; IDL port of numpy's unique() that returns the inverse
;;; i.e., if _, result = np.unique(arr, return_inverse=True)
;;; then result is the same as result = invert_uniq(arr)
;;; Written by Daniel J. Stevens (UM_Duluth), 2024 July 23
	idx = uniq(arr, sort(arr))
	idx_inv = uintarr(n_elements(arr))
	for i=0, n_elements(arr[idx])-1 do begin
		tmp_idx = where(arr eq arr[idx[i]])
		idx_inv[tmp_idx] = i
	endfor
	return, idx_inv
 end