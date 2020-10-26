pro fithat3notran, debug=debug, verbose=verbose, maxsteps=maxsteps, nthin=nthin

if n_elements(nthin) eq 0 then nthin=1
path = filepath('',root_dir=getenv('EXOFAST_PATH'),subdir=['examples','hat3'])

;; Fit using the Torres relation
exofastv2, nplanets=1, $;tranpath=path+'n20070428.Sloani.KepCam.dat', $
           rvpath=path+'HAT-3b.*.rv',priorfile=path+'hat3.torres.priors', $
           prefix=path+'fitresults' + path_sep()+'HAT-3b.notran.', maxsteps=maxsteps, $
           nthin=nthin, circular=[1], fitrv=[1],$;fittran=[1], $
           debug=debug, verbose=verbose, /nomist, /torres

end
