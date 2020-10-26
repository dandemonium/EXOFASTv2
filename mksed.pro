;+
; NAME:
;   STROM_CONV
; PURPOSE:
;    Translates stromgren color combinations from the catalog to
;    individual uvby magnitudes
; Modification 
;    2018-04-12: Jason Eastman, CfA
;                Renamed, documented, and cleaned up for distribution with EXOFASTv2
;
;-
function strom_conv,V,sigV,by,sigby,m1,sigm1,c1,sigc1,silent=silent

if n_params() lt 8 then begin
  print,'syntax: result=strom_conv(V,sigV,by,sigby,m1,sigm1,c1,sigc1)'
  retall
endif

y = V
b = V + by 
u = V + 3*by + 2*m1 + c1
v = V + 2*by + m1

sigy = sigV
sigb = sqrt(sigV^2 + sigby^2)
sigu = sqrt(sigV^2 + (3*sigby)^2 + (2*sigm1)^2 + sigc1^2)
sigv = sqrt(sigV^2 + (2*sigby)^2 + sigm1^2)

if not keyword_set(silent) then begin
  print,'uvby:',u,v,b,y
  print,'sig_uvby:',sigu,sigv,sigb,sigy
endif

return,[u,sigu,v,sigv,b,sigb,y,sigy]

end

;+
; NAME:
;   CROSSREF
; PURPOSE:
;     Finds the closest match to refra and refdec from an array of RAs
;     and Decs.
; Modification 
;    2018-04-12: Jason Eastman, CfA
;                Renamed, documented, and cleaned up for distribution with EXOFASTv2
;
;-
function crossref, refra, refdec, raarr, decarr

sep = angsep(refra*!dpi/180d0+dblarr(n_elements(raarr)),$
             refdec*!dpi/180d0+dblarr(n_elements(raarr)),$
             raarr*!dpi/180d0,decarr*!dpi/180d0)
junk = min(sep,match)
return, match

end

;+
; NAME:
;    MKSED
; PURPOSE:
;    Given an object name (resolved with simbad) or RA and DEC, return
;    available broad band photometry from a variety of catalogs in the
;    format expected by EXOFASTv2 to perform an SED fit.
;    $EXOFAST_PATH/sed/mag2fluxconv.pro reads this file.
;
;    NOTE: Priors on extinction and parallax should be included (see
;    mkpriors.pro), as well as a prior on Teff.
;
; INPUTS:
; 
;    References included in the output file -- please cite the appropriate catalogs.
; Modification 
;    20??-??-??: Written: Keivan Stassun, Vanderbilt
;    2018-04-12: Jason Eastman, CfA
;                Renamed, documented, and cleaned up for distribution
;                with EXOFASTv2
;    2019-01-18: Now can be run without a license (use mksed.sav)
;
;-
pro mksed,name,outfile,dist=dist,galdist=galdist,$
          apass=apass,merm=merm,stromgren=stromgren,kepler=kepler,galex=galex,$
          nogaia=nogaia,nowise=nowise,over=over,cfa=cfa,$
          ra=ra,dec=dec,silent=silent, parallax=parallax, uparallax=uparallax

;; for use without a license
if lmgr(/vm) or lmgr(/runtime) then begin
   par = command_line_args(count=numargs)
   
   name = par[0]
   if numargs gt 1 then outfile = par[1]

   ;; set keywords
   if (where(strupcase(par) eq 'APASS'))[0] ne -1 then apass = 1
   if (where(strupcase(par) eq 'MERM'))[0] ne -1 then merm = 1
   if (where(strupcase(par) eq 'STROMGREN'))[0] ne -1 then stromgren = 1
   if (where(strupcase(par) eq 'KEPLER'))[0] ne -1 then kepler = 1
   if (where(strupcase(par) eq 'GALEX'))[0] ne -1 then galex = 1
   if (where(strupcase(par) eq 'NOGAIA'))[0] ne -1 then nogaia = 1
   if (where(strupcase(par) eq 'NOWISE'))[0] ne -1 then nowise = 1
   if (where(strupcase(par) eq 'OVER'))[0] ne -1 then over = 1
   if (where(strupcase(par) eq 'CFA'))[0] ne -1 then cfa = 1
   if (where(strupcase(par) eq 'SILENT'))[0] ne -1 then silent = 1

   for i=0L, numargs-1 do begin
      if strpos(par[i],'=') ne -1 then begin
         entries = strsplit(par[i],'=',/extract)
         if strupcase(entries[0]) eq 'DIST' then dist = double(entries[1])
         if strupcase(entries[0]) eq 'GALDIST' then galdist = double(entries[1])
         if strupcase(entries[0]) eq 'RA' then ra = double(entries[1])
         if strupcase(entries[0]) eq 'DEC' then dec = double(entries[1])
         if strupcase(entries[0]) eq 'PARALLAX' then parallax = double(entries[1])
         if strupcase(entries[0]) eq 'UPARALLAX' then uparallax = double(entries[1])
      endif
   endfor

endif

if n_elements(parallax) eq 0 then parallax = 0d0
if n_elements(uparallax) eq 0 then uparallax = -1d0

star=name
if n_elements(dist) eq 0 then dist=2d0
if n_elements(galdist) eq 0 then galdist=dist
if n_elements(ra) eq 0 or n_elements(dec) eq 0 then querysimbad,star,ra,dec,cfa=cfa else star=[ra,dec]

if n_elements(ra) eq 0 or n_elements(dec) eq 0 then begin
   message, 'Star name not resolved by SIMBAD. Sometimes their server is down and trying again helps. Otherwise, you must supply the RA and Dec directly.',/continue
   return
endif

;; if outfile not specified, print results to the screen
if n_elements(outfile) eq 0 then lun = -1 $
else begin
   if not keyword_set(over) and file_test(outfile) then begin
      print,'File '+outfile+' exists and over keyword not specified. Quitting.'
      return
   endif
   openw,lun,outfile, /get_lun
endelse

printf, lun, '# Band, Mag, used mag error, catalog mag error'
fmt = '(a-6,x,f7.3,x,f5.3,x,f5.3)'

;; get GALEX -- by default, skip; models are unreliable here
qgalex=Exofast_Queryvizier('II/312/ais',star,galdist/60d0,/silent,/all,cfa=cfa)
if long(tag_exist(qgalex,'fuv',/quiet)) ne 0L then begin
   printf, lun, '# Galex DR5, Bianchi+ 2011'
   printf, lun, '# http://adsabs.harvard.edu/abs/2011Ap%26SS.335..161B'
   if n_elements(qgalex) gt 1 then begin
      print,'Warning: More than 1 GALEX source found; using nearest one only.'
      printf, lun,'# Warning: More than 1 GALEX source found; using nearest one only.'
      junk = min(qgalex._r,m)
      qgalex=qgalex[m]
   endif
   if (qgalex.fuv gt -99 and finite(qgalex.e_fuv)) or (qgalex.nuv gt -99 and finite(qgalex.e_nuv)) then begin
      printf, lun,'# including the GALEX is likely to bias the result due to model errors. Use with caution.'
   endif
   if qgalex.fuv gt -99 and finite(qgalex.e_fuv) then begin
      line = string('galFUV',qgalex.fuv,max([0.1d,qgalex.e_fuv]),qgalex.e_fuv, format=fmt) 
      if keyword_set(galex) then printf, lun, line $
      else printf, lun, '#' + line
   endif
   if qgalex.nuv gt -99 and finite(qgalex.e_nuv) then begin
      line = string('galNUV',qgalex.nuv,max([0.1d,qgalex.e_nuv]),qgalex.e_nuv, format=fmt)
      if keyword_set(galex) then printf, lun,line $
      else printf, lun, '#' + line
   endif
endif else qgalex={fuv_6:-99.,nuv_6:-99.}


; get Mermilliod 1991 UBV ;; by default skip (Gaia better)
qmermilliod91=Exofast_Queryvizier('II/168/ubvmeans',star,dist/60.,/silent,/all,cfa=cfa) 
if long(tag_exist(qmermilliod91,'vmag',/quiet)) ne 0L then begin
   if n_elements(qmermilliod91) gt 1 then begin
      print,'Warning: More than 1 Mermilliod source found; using nearest one only.'
      printf, lun,'# Warning: More than 1 Mermilliod source found; using nearest one only.'
      junk = min(qmermilliod91._r,m)
      qmermilliod91=qmermilliod91[m]
   endif
   
   printf, lun, '# Mermilliod, 1994'
   printf, lun, '# http://adsabs.harvard.edu/abs/1994yCat.2193....0M'
   if qmermilliod91.u_b+qmermilliod91.b_v+qmermilliod91.vmag gt -9 and $
      finite(qmermilliod91.e_u_b) and finite(qmermilliod91.e_b_v) and $
      finite(qmermilliod91.e_vmag) then begin
      line = string('U',qmermilliod91.u_b+qmermilliod91.b_v+qmermilliod91.vmag,max([0.02d,sqrt(qmermilliod91.e_u_b^2+qmermilliod91.e_b_v^2+qmermilliod91.e_vmag^2)]), format=fmt)
      if keyword_set(merm) then printf, lun, line $
      else printf, lun, '#' + line
   endif
   
   if qmermilliod91.b_v+qmermilliod91.vmag gt -9 and $
      finite(qmermilliod91.e_b_v) and $
      finite(qmermilliod91.e_vmag) then begin
      line = string('B',qmermilliod91.b_v+qmermilliod91.vmag,max([0.02d,sqrt(qmermilliod91.e_b_v^2+qmermilliod91.e_vmag^2)]), format=fmt)
      if keyword_set(merm) then printf, lun, line $
      else printf, lun, '#' + line
   endif
   if qmermilliod91.vmag gt -9 and finite(qmermilliod91.e_vmag) then begin
      line = string('V',qmermilliod91.vmag,max([0.02d,qmermilliod91.e_vmag]), format=fmt)
      if keyword_set(merm) then printf, lun, line $
      else printf, lun, '#' + line
   endif
endif else qmermilliod91={vmag:-99.}

; Tycho-2
qtyc2=Exofast_Queryvizier('I/259/TYC2',star,dist/60.,/silent,/all,cfa=cfa)
if long(tag_exist(qtyc2,'BTMAG',/quiet)) ne 0L then begin
   printf, lun, '# Tycho catalog, Hoeg+ 2000'
   printf, lun, '# http://adsabs.harvard.edu/abs/2000A%26A...355L..27H'
   if n_elements(qtyc2) gt 1 then begin
      print,'Warning: More than 1 Tycho-2 source found; using nearest one only.'
      printf, lun,'# Warning: More than 1 Tycho-2 source found; using nearest one only.'
      junk = min(qtyc2._r,m)
      qtyc2=qtyc2[m]
   endif
   if qtyc2.btmag gt -9 and finite(qtyc2.e_btmag) then printf, lun,'BT',qtyc2.btmag,max([0.02d,qtyc2.e_btmag]),qtyc2.e_btmag, format=fmt
   if qtyc2.vtmag gt -9 and finite(qtyc2.e_vtmag) then printf, lun,'VT',qtyc2.vtmag,max([0.02d,qtyc2.e_vtmag]),qtyc2.e_vtmag, format=fmt
endif else qtyc2={btmag:-99.,vtmag:-99.}

; Paunzen, 2015 Stromgren
qpaunzen15=Exofast_Queryvizier('J/A+A/580/A23/catalog',star,dist/60.,/silent,/all,cfa=cfa)
if long(tag_exist(qpaunzen15,'vmag',/quiet)) ne 0L then begin
   if n_elements(qpaunzen15) gt 1 then begin
      print,'Warning: More than 1 Paunzen source found; using nearest one only.'
      printf, lun,'# Warning: More than 1 Paunzen source found; using nearest one only.'
      junk = min(qpaunzen15._r,m)
      qpaunzen15=qpaunzen15[m]
   endif
   ubvymags=strom_conv(qpaunzen15.vmag,max([0.01d,qpaunzen15.e_vmag]),$
                       qpaunzen15.b_y,max([0.02d,qpaunzen15.e_b_y]),$
                       qpaunzen15.m1,max([0.02d,qpaunzen15.e_m1]),$
                       qpaunzen15.c1,max([0.02d,qpaunzen15.e_c1]),/silent)
   printf, lun, '# Stromgren photometry, Paunzen, 2015'
   printf, lun, '# http://adsabs.harvard.edu/abs/2015A%26A...580A..23P'
   if keyword_set(stromgren) then comment = '' else comment = '#'

   if ubvymags(0) gt -9 then printf, lun,comment+'uStr',ubvymags(0),max([0.02d,ubvymags(1)]), format=fmt
   if ubvymags(2) gt -9 then printf, lun,comment+'vStr',ubvymags(2),max([0.02d,ubvymags(3)]), format=fmt
   if ubvymags(4) gt -9 then printf, lun,comment+'bStr',ubvymags(4),max([0.02d,ubvymags(5)]), format=fmt
   if ubvymags(6) gt -9 then printf, lun,comment+'yStr',ubvymags(6),max([0.02d,ubvymags(7)]), format=fmt
endif else qpaunzen15={vmag:-99.}

; UCAC4 - check if APASS BV is actually just TYC2 BV, and adjust errors from integer to centimag
qucac4=Exofast_Queryvizier('UCAC4',star,dist/60.,/silent,/all,cfa=cfa)
if long(tag_exist(qucac4,'bmag',/quiet)) ne 0L then begin
   printf, lun,'# UCAC4, Zacharias+, 2012'
   printf, lun,'# http://adsabs.harvard.edu/abs/2012yCat.1322....0Z'
   if n_elements(qucac4) gt 1 then begin
      print,'Warning: More than 1 UCAC-4 source found; using nearest one only.'
      printf, lun,'# Warning: More than 1 UCAC-4 source found; using nearest one only.'
      junk = min(qucac4._r,m)
      qucac4=qucac4[m]
   endif
   if keyword_set(apass) then begin
      printf, lun, '# APASS DR6 (via UCAC4), Henden+ 2016'
      printf, lun, '# http://adsabs.harvard.edu/abs/2016yCat.2336....0H'
      if qucac4.bmag ne qtyc2.btmag and qucac4.bmag gt -9 and qucac4.e_bmag ne 99 then printf, lun,'B',qucac4.bmag,max([0.02d,qucac4.e_bmag*0.01d]),qucac4.e_bmag*0.01d, format=fmt
      if qucac4.vmag ne qtyc2.vtmag and qucac4.bmag gt -9 and qucac4.e_vmag ne 99 then printf, lun,'V',qucac4.vmag,max([0.02d,qucac4.e_vmag*0.01d]),qucac4.e_vmag*0.01d, format=fmt
      if qucac4.gmag gt -9 then printf, lun,'gSDSS',qucac4.gmag,max([0.02d,qucac4.e_gmag*0.01d]),qucac4.e_gmag*0.01d, format=fmt
      if qucac4.rmag gt -9 then printf, lun,'rSDSS',qucac4.rmag,max([0.02d,qucac4.e_rmag*0.01d]),qucac4.e_rmag*0.01d, format=fmt
      if qucac4.imag gt -9 then printf, lun,'iSDSS',qucac4.imag,max([0.02d,qucac4.e_imag*0.01d]),qucac4.e_imag*0.01d, format=fmt
   endif
   if qucac4.Jmag gt -9  or qucac4.Hmag gt -9 or qucac4.Kmag gt -9 then begin
      printf, lun, '# 2MASS (via UCAC4), Cutri+ 2003'
      printf, lun, '# http://adsabs.harvard.edu/abs/2003yCat.2246....0C'   
   endif
   if qucac4.Jmag gt -9 then printf, lun,'J2M',qucac4.Jmag,max([0.02d,qucac4.e_Jmag]),qucac4.e_Jmag, format=fmt
   if qucac4.Hmag gt -9 then printf, lun,'H2M',qucac4.Hmag,max([0.02d,qucac4.e_Hmag]),qucac4.e_Hmag, format=fmt
   if qucac4.Kmag gt -9 then printf, lun,'K2M',qucac4.Kmag,max([0.02d,qucac4.e_Kmag]),qucac4.e_Kmag, format=fmt
endif else qucac4={Jmag:-99.,Hmag:-99.,Kmag:-99.}

;; ALLWISE
qwise=Exofast_Queryvizier('II/328/allwise',star,dist/60.,/silent,cfa=cfa)
if long(tag_exist(qwise,'w1mag',/quiet)) ne 0L then begin
   printf, lun, '# AllWISE, Cutri+, 2013'
   printf, lun, '# http://adsabs.harvard.edu/abs/2013yCat.2328....0C'
   if n_elements(qwise) gt 1 then begin
      print,'Warning: More than 1 WISE source found; using nearest one only.'
      printf, lun,'# Warning: More than 1 WISE source found; using nearest one only.'
      match = crossref(ra,dec,qwise.raj2000,qwise.dej2000)
      qwise = qwise[match]
   endif
   if keyword_set(nowise) then comment = '#' else comment = ''

   if qwise.w1mag gt -9 and finite(qwise.e_w1mag) then printf, lun,comment + 'WISE1',qwise.w1mag,max([0.03d,qwise.e_w1mag]),qwise.e_w1mag, format=fmt
   if qwise.w2mag gt -9 and finite(qwise.e_w2mag) then printf, lun,comment + 'WISE2',qwise.w2mag,max([0.03d,qwise.e_w2mag]),qwise.e_w2mag, format=fmt
   if qwise.w3mag gt -9 and finite(qwise.e_w3mag) then printf, lun,comment + 'WISE3',qwise.w3mag,max([0.03d,qwise.e_w3mag]),qwise.e_w3mag, format=fmt
   if qwise.w4mag gt -9 and finite(qwise.e_w4mag) then printf, lun,comment + 'WISE4',qwise.w4mag,max([0.1d,qwise.e_w4mag]),qwise.e_w4mag, format=fmt

   if long(tag_exist(qucac4,'Jmag',/quiet)) eq 0L then begin ; this means that we didn't get 2MASS via UCAC so get them from WISE
      ;; does this ever happen?
      printf, lun, '# 2MASS (via WISE), Cutri+ 2003'
      printf, lun, '# http://adsabs.harvard.edu/abs/2003yCat.2246....0C'
      if qwise.Jmag gt -9 then printf, lun,'J2M',qwise.Jmag,max([0.02d,qwise.e_Jmag]),qwise.e_Jmag, format=fmt
      if qwise.Hmag gt -9 then printf, lun,'H2M',qwise.Hmag,max([0.02d,qwise.e_Hmag]),qwise.e_Hmag, format=fmt
      if qwise.Kmag gt -9 then printf, lun,'K2M',qwise.Kmag,max([0.02d,qwise.e_Kmag]),qwise.e_Kmag, format=fmt
   endif else if qucac4.jmag lt -9 or qucac4.hmag lt -9 or qucac4.kmag lt -9 then begin
      printf, lun, '# 2MASS (via WISE), Cutri+ 2003'
      printf, lun, '# http://adsabs.harvard.edu/abs/2003yCat.2246....0C'
      if qucac4.jmag lt -9 and qwise.Jmag gt -9 then printf, lun,'J2M',qwise.Jmag,max([0.02d,qwise.e_Jmag]),qwise.e_Jmag, format=fmt
      if qucac4.hmag lt -9 and qwise.Hmag gt -9 then printf, lun,'H2M',qwise.Hmag,max([0.02d,qwise.e_Hmag]),qwise.e_Hmag, format=fmt
      if qucac4.kmag lt -9 and qwise.Kmag gt -9 then printf, lun,'K2M',qwise.Kmag,max([0.02d,qwise.e_Kmag]),qwise.e_Kmag, format=fmt
   endif

endif else qwise={Jmag:-99.,Hmag:-99.,Kmag:-99.}

; 2MASS (if not from UCAC or WISE)
if (qucac4.jmag lt -90 or ~finite(qucac4.jmag)) and (qwise.Jmag lt -90 or ~finite(qwise.Jmag)) then begin
   printf, lun, '# 2MASS, Cutri+ 2003'
   printf, lun, '# http://adsabs.harvard.edu/abs/2003yCat.2246....0C'
   q2mass=Exofast_Queryvizier('II/246/out',star,dist/60.,/silent,cfa=cfa)
   if long(tag_exist(q2mass,'Jmag',/quiet)) ne 0L then begin
      if n_elements(q2mass) gt 1 then begin
         print,'Warning: More than 1 2mass source found; using nearest one only.'
         printf, lun,'# Warning: More than 1 2mass source found; using nearest one only.'
         match = crossref(ra,dec,q2mass.raj2000,q2mass.dej2000)
         q2mass = q2mass[match]
      endif
      if q2mass.Jmag gt -9 then printf, lun,'J2M',q2mass.Jmag,max([0.02d,q2mass.e_Jmag]),q2mass.e_Jmag, format=fmt
      if q2mass.Hmag gt -9 then printf, lun,'H2M',q2mass.Hmag,max([0.02d,q2mass.e_Hmag]),q2mass.e_Hmag, format=fmt
      if q2mass.Kmag gt -9 then printf, lun,'K2M',q2mass.Kmag,max([0.02d,q2mass.e_Kmag]),q2mass.e_Kmag, format=fmt
   endif
endif

;; Gaia DR2
qgaia=Exofast_Queryvizier('I/345/gaia2',star,dist/60.,/silent,cfa=cfa,/all)
if long(tag_exist(qgaia,'gmag',/quiet)) ne 0L then begin
   printf, lun, '# Gaia DR2, Gaia Collaboration, 2018'
   printf, lun, '# http://adsabs.harvard.edu/abs/2016A%26A...595A...2G'
   if n_elements(qgaia) gt 1 then begin
      print,'Warning: More than 1 Gaia source found; using nearest one only.'
      printf, lun,'# Warning: More than 1 Gaia source found; using nearest one only.'
      match = crossref(ra, dec, qgaia.ra_icrs,qgaia.de_icrs)
      qgaia = qgaia[match]
   endif
   
;      print, name
;      print, 'parallax',qgaia.plx, qgaia.e_plx, ' # Gaia DR2'
   parallax = qgaia.plx
   uparallax = qgaia.e_plx
   
   if keyword_set(nogaia) then comment = '#' else comment = ''
   
   if qgaia.gmag gt -9 and finite(qgaia.e_gmag) then printf, lun,comment+'Gaia',qgaia.gmag,max([0.02d,qgaia.e_gmag]),qgaia.e_gmag, format=fmt
   if qgaia.bpmag gt -9 and finite(qgaia.e_bpmag) then printf, lun,comment+'GaiaBP',qgaia.bpmag,max([0.02d,qgaia.e_bpmag]),qgaia.e_bpmag, format=fmt
   if qgaia.rpmag gt -9 and finite(qgaia.e_rpmag) then printf, lun,comment+'GaiaRP',qgaia.rpmag,max([0.02d,qgaia.e_rpmag]),qgaia.e_rpmag, format=fmt
   
endif

; KIS DR2
qkis=Exofast_Queryvizier('J/AJ/144/24/kisdr2',star,dist/60.,/silent,cfa=cfa,/all)
if long(tag_exist(qkis,'KIS',/quiet)) ne 0L then begin
   printf, lun, '# KIS DR2, Greiss+ 2012'
   printf, lun, '# http://adsabs.harvard.edu/abs/2012AJ....144...24G'
   if n_elements(qkis) gt 1 then begin
      print,'Warning: More than 1 KIS source found; using nearest one only.'
      printf, lun,'# Warning: More than 1 KIS source found; using nearest one only.'
      match = crossref(ra, dec, qkis.raj2000, qkis.dej2000)
      qkis = qkis[match]
   endif
   
   if keyword_set(kepler) then comment = '' else comment = '#'
   
   if qkis.umag gt -9 and finite(qkis.e_umag) then printf, lun,comment+'UKIS',qkis.umag,max([0.02d,qkis.e_umag]),qkis.e_umag, format=fmt
   if qkis.gmag gt -9 and finite(qkis.e_gmag) then printf, lun,comment+'gKIS',qkis.gmag,max([0.02d,qkis.e_gmag]),qkis.e_gmag, format=fmt
   if qkis.rmag gt -9 and finite(qkis.e_rmag) then printf, lun,comment+'rKIS',qkis.rmag,max([0.02d,qkis.e_rmag]),qkis.e_rmag, format=fmt
   if qkis.imag gt -9 and finite(qkis.e_imag) then printf, lun,comment+'iKIS',qkis.imag,max([0.02d,qkis.e_imag]),qkis.e_imag, format=fmt
endif else qkis={umag:-99.,gmag:-99.,rmag:-99.,imag:-99.}

if lun ne -1 then free_lun, lun

end

