pro mkgaiaspecphot,ticid=ticid,dist=dist,gaiaspfile=gaiaspfile

if n_elements(gaiaspfile) eq 0 then gaiaspfile = ticid + '.Gaia.sed'


qgaiasp=Exofast_Queryvizier('I/355/xpsample',ticid,dist/60.,/silent,cfa=cfa,/all)
if n_elements(qgaiasp) gt 1 then begin ; multiple matches, so take closest
  index = where(qgaiasp[*]._r eq min(qgaiasp[*]._r))
  qgaiasp = qgaiasp[index]
endif
;; Gaia lambda in nm, Gaia flux in W/m^2/Hz
;; SED plot in microns, erg/s/cm^2
;; interpolate onto atmosphere wavelength scale now?
flux = qgaiasp.flux*1d3*qgaiasp.lambda ;; W/m^2/Hz -> erg/s/cm^2
fluxerr = qgaiasp.e_flux*1d3*qgaiasp.lambda ;; W/m^2/Hz -> erg/s/cm^2
lambda = qgaiasp.lambda/1d3 ;; nm -> um
exofast_forprint, lambda, flux, fluxerr, textout=gaiaspfile,/nocomment,format='(f5.3,x,e12.6,x,e12.6)'

print, 'Successfully retrieved Gaia DR3 spectrophotometry information. See ' + gaiaspfile

end