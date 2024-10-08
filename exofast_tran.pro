;+
; NAME:
;   EXOFAST_TRAN
;
;   PURPOSE:

;      Computes a transit model given a complete set of physical
;      parameters
;
;  CALLING SEQUENCE:
;      model = EXOFAST_TRAN(time, inc, ar, tp, period, e, omega, p, u1, u2, f0, $
;                       rstar=rstar, thermal=thermal, reflect=reflect, $
;                       dilute=dilute, tc=tc, q=q,x1=x1,y1=y1,z1=z1,
;                       au=au,c=c)
;
;  INPUTS:
;     TIME   - The BJD_TDB time of the model to compute. Scalar or
;              array.
;     INC    - The inclination of the planetary orbit, in radians
;     AR     - a/Rstar, the semi-major axis of the planetary orbit in
;              units of stellar radii
;     TP     - The time of periastron, in BJD_TDB
;     PERIOD - The planetary orbital period, in days
;     E      - The planetary eccentricity
;     OMEGA  - The argument of periastron, in radians
;     P      - Rp/Rstar, the planetary radius in units of stellar
;              radii
;     U1     - The linear limb darkening parameter
;     U2     - The quadratic limb darkening parameter
;     F0     - The transit baseline level
;
; OPTIONAL INPUTS:
;     RSTAR   - The stellar radius **in AU**. If supplied, it will
;               apply the light travel time correction to the target's
;               barycenter
;     THERMAL - The thermal emission contribution from the planet, in ppm.
;     REFLECT - The reflected light from the planet, in ppm. 
;     DILUTE  - The fraction of to basline flux that is due to
;               contaminating sources. DILUTE = F2/(F1+F2), where F1
;               is the flux from the host star, and F2 is the flux
;               from all other sources in the aperture.
;ELLIPSOIDAL  - The flux variations caused by tidal ellipsoidal variations, in ppm.
;     TC      - The time of conjunction, used for calculating the
;               phase of the reflected light component. If not
;               supplied, it is computed from e, omega, period.
;     Q       - M1/M2. The mass ratio of the primary to companion. If
;               supplied, the stellar reflex motion is computed and
;               applied.
;     X1      - The X motion of the star due to all other bodies in the
;               system, in units of Rstar with the origin at the
;               barycenter. Output from exofast_getb2.pro.
;     Y1      - Same as X1, but in the Y direction
;     Z1      - Same as X1, but in the Z direction
;     AU      - The value of the AU, in solar radii, default is
;               215.094177d0
;     C       - The speed of light, in AU/day (default is computed in
;               TARGET2BJD).
;  OUTPUTS:
;    MODEL - The transit model as a function of time
; 
;  REVISION HISTORY:
;    2015 (?) - Written by Jason Eastman (CfA)
;    2018/10  - Documented (JDE)
;    2019/01/28 - Replaced incorrect arg_present check with n_elements
;               check. Didn't convert to target frame before,
;               as called by exofast_chi2v2.pro
;    2023/03/16 - Adds eBEER formulae instead of DJS formulae
;-
function exofast_tran, time, inc, ar, tp, period, e, omega, p, u1, u2, f0, $
                       rstar=rstar, thermal=thermal, reflect=reflect, beam=beam, ellipsoidal=ellipsoidal, $
                       phaseshift=phaseshift, $
                       dilute=dilute, tc=tc, q=q, au=au,c=c;,x1=x1,y1=y1,z1=z1

if n_elements(thermal) eq 0 then thermal = 0d0
if n_elements(reflect) eq 0 then reflect = 0d0
if n_elements(dilute) eq 0 then dilute = 0d0
if n_elements(ellipsoidal) eq 0 then ellipsoidal = 0d0
if n_elements(beam) eq 0 then beam = 0d0
if n_elements(AU) eq 0 then AU = 215.094177d0
if n_elements(phaseshift) eq 0 then phaseshift = 0d0

;; if we have the stellar radius, we can convert time to the
;; target's barycentric frame
if n_elements(rstar) ne 0 then begin
   transitbjd = bjd2target(time, inclination=inc, a=ar*rstar, tp=tp, $
                           period=period, e=e, omega=omega,q=q,c=c)
endif else transitbjd = time

;; the impact parameter for each BJD
z = exofast_getb2(transitbjd, i=inc, a=ar, tperiastron=tp, period=period,$
                     e=e,omega=omega,z2=depth,x2=x,y2=y,q=q)

ntime = n_elements(time)


if reflect ne 0d0 or ellipsoidal ne 0d0 or beam ne 0d0 then begin 
   if e ne 0 then begin
      meananom = 2.d0*!dpi*(1.d0 + (transitbjd - tp)/period mod 1)
      ;; if eccentricitys given, integrate the orbit
      if n_elements(e) ne 0 then begin
         if n_elements(omega) eq 0 then message, $
            'ERROR: omega must be specified if e is specified'
         if e lt 0 then begin
            e = -e
            omega = omega + !dpi
         endif else begin
            e = e
            omega = omega
         endelse
	 if omega < 0 then omega += !dpi/2d0
         eccanom = exofast_keplereq(meananom, e)
         trueanom = 2.d0*atan(sqrt((1.d0 + e)/(1.d0 - e))*tan(eccanom/2.d0))
      endif else begin
         e=0.d0
      endelse
   ;; standard definition of omega for circular orbits
   endif else begin
     e = e
	trueanom = 0d0	
	if n_elements(omega) eq 0 then omega = !dpi/2.d0
   endelse
   sep = (1.-(e*e))/(1.+(e*cos(trueanom)));*q/(1.+q)
   beta_ebeer = 1./sep
endif


;; Primary transit
modelflux = dblarr(n_elements(time))+1d0
primary = where(depth lt 0, complement=secondary)
if primary[0] ne - 1 then begin
   exofast_occultquad_cel, z[primary], u1, u2, p, mu1
   modelflux[primary] =  mu1
endif

;; calculate the fraction of the planet that is visible for each time
if thermal ne 0d0 or reflect ne 0d0 then begin
   planetvisible = dblarr(n_elements(time)) + 1d0
   if secondary[0] ne - 1 then begin
      exofast_occultquad_cel, z[secondary]/abs(p), 0, 0, 1d0/p, mu1
      planetvisible[secondary] = mu1
   endif
endif

;; thermal emission from planet (isotropic)
if thermal ne 0d0 then modelflux += 1d-6*thermal*planetvisible

;if hotspot ne 0d0 then begin
;   ;; This makes flux=0 during primary transit (Thanks Sam Quinn!)
;   modelflux-=(1d-6*hotspot/2d0)*(cos(2d0*!dpi*(transitbjd-tc0)/period+phaseshift/360d0)-1d0)*planetvisible
;endif   

;; phase-dependent reflection off planet
if reflect ne 0d0 then begin

   ;; This is non-physical! negative flux during primary transit
   ;modelflux-=1d-6*reflect*cos(2d0*!dpi*(transitbjd-tc0)/period)*planetvisible
  if e eq 0 then begin  ;;; commented out June 6, 2023;
    if n_elements(tc) eq 0 then begin
      phase = exofast_getphase(e,omega,/primary)  
      tc0 = tp - phase*period
    endif else tc0 = tc
   ;; This makes flux=0 during primary transit (Thanks Sam Quinn!)
   modelflux-=(1d-6*reflect/2d0)*(cos(2d0*!dpi*(transitbjd-tc)/period)-1d0)*planetvisible
  endif else begin
   ;;modelflux+=1d-6*reflect*sin(!dpi*(transitbjd-tc0)/period)^2d0*planetvisible
;  if e eq 0 then begin  ;;; commented out June 6, 2023
;    cosz = -1d0*sin(inc)*cos(2d0*!dpi*(transitbjd-tc)/period)
;    lambertamp = reflect*(sin(acos(cosz)) + (!dpi - acos(cosz))*cosz)/!dpi
;    modelflux += 1d-6*lambertamp*planetvisible ;; assumes the sin(i)
 ;;;terms in A_ref  propto sin(i) and Psi(lam)/sin(i) are cancelled out already.
;  endif else begin
   modelflux += 1d-6*reflect*(beta_ebeer^(2d0))*(0.64 - (sin(inc)*sin(omega+trueanom))+(0.18*sin(inc)*sin(inc)*(1.-cos(2.*(omega+trueanom)))))*planetvisible ;; this is the eBEER formalism, but with a * that I omitted.
   endelse
endif

;; add beaming and ellipsoidal variations
if ellipsoidal ne 0d0 then begin
  if e eq 0 then begin  ;;; commented out June 6, 2023;
    if n_elements(tc) eq 0 and n_elements(tc0) eq 0 then begin
      phase = exofast_getphase(e,omega,/primary)  
      tc0 = tp - phase*period
    endif else tc0 = tc
  modelflux -= ellipsoidal*1d-6*cos(2d0*!dpi*(transitbjd-tc)/(period/2d0))
  endif else begin
    modelflux+= 1d-6*ellipsoidal*cos(2.*(omega+trueanom))*beta_ebeer^(3.)
  endelse
endif

if beam ne 0d0 then begin
   if n_elements(tc) eq 0 and n_elements(tc0) eq 0 then begin
      phase = exofast_getphase(e,omega,/primary)  
      tc0 = tp - phase*period
   endif else tc0 = tc
   if e eq 0 then modelflux += beam*1d-6 * sin(2d0*!dpi*(transitbjd-tc)/period) else $
   modelflux += -1d-6*beam*sin(inc)*cos(trueanom+omega)/sqrt(1.-(e*e))
endif

;; normalization and dilution due to neighboring star
if dilute ne 0d0 then modelflux = f0*(modelflux*(1d0-dilute)+dilute) $
else modelflux *= f0


return, modelflux

end
