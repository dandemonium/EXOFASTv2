;+
; NAME:
;   exofast_plotdist_corner
; PURPOSE:
;   Plots the distributions from MCMC chains (Probability Distribution
;   Functions (PDFs) and covariances; returns median values and 68%
;   confidence intervals.
;
; DESCRIPTION:
;
; CALLING SEQUENCE:
;   exofast_plotdist_corner, pars [,MEDIANPARS, PARNAMES=, UNITS=,
;                     /NOCOVAR, PDFNAME=, COVARNAME=, BESTPARS=, PROBS=, /DEGREES]
;
; INPUTS:
;   SS         - The stellar structure created by EXOFASTv2, populated
;                with an array of values for each parameter
;
; OPTIONAL INPUTS:
;   PDFNAME    - The filename of the output postscript plot of
;                PDFs. Default is pdf.ps
;   COVARNAME  - The filename of the output postscript plot of
;                covariances of fitted parameters. Default is covar.ps. 
;   PROBS      - Array of probability contours to draw on the covariance
;                plots. Default is erf([1d,2d]/sqrt(2d0)) ~ [0.68,0.95].
;   LOGNAME    - The name of the logfile to print messages to.
;
; OPTIONAL KEYWORDS:
;   NOCOVAR    - If set, no covariances plots will be made. Since
;                there will be NFIT*(NFIT-1)/2 plots, this can be
;                time consuming.
; OUTPUTS:
;   MEDIANPARS - A 3 x NPARS array of median values, +err, -err (68%
;                confidence interval)
;   PDFNAME    - A postscript plot of Probablity Distribution
;                Functions (PDFs)
;   COVARNAME  - A postscript plot of covariances.
;
; REVISION HISTORY:
;   2009/11/24 - Written: Jason Eastman - The Ohio State University
;   2013/03/12 - Fixed min/max, made slightly more efficient
;   2018/02/13 - Covariance plot now a corner plot inspired by
;                corner.py (much prettier)
;   2018/03/15 - Identify and discard bad chains before inference
;   2018/10/18 - Remove misleading (and unused) ANGULAR parameter;
;                functionality applied with unit property in
;                summarizepar.
;-
pro exofast_plotdist_corner, ss, nocovar=nocovar, $
                             pdfname=pdfname, covarname=covarname, $
                             probs=probs,logname=logname, $
                             useallpars=useallpars, csvfile=csvfile, mask=mask

if n_elements(pdfname) eq 0 then pdfname = 'pdf.ps'
if n_elements(covarname) eq 0 then covarname = 'covar.ps'

;; 68% and 95% probability contours
if n_elements(probs) eq 0 then probs = erf([1d,2d]/sqrt(2d0))

nsteps = ss.nsteps/ss.nchains
chi2 = reform((*ss.chi2),nsteps,ss.nchains)
burnndx = ss.burnndx
goodchains=(*ss.goodchains)
ngoodchains = n_elements(goodchains)

;; prepare the postscript device
mydevice=!d.name
set_plot, 'PS'
aspect_ratio=1.5
xsize = 18
ysize=xsize/aspect_ratio
!p.font=0

defsysv, '!GDL', exists=runninggdl  
if ~runninggdl then device, filename=pdfname, set_font='Times',/tt_font, xsize=xsize,ysize=ysize,/color,bits=24
charsize=1.5
loadct,39,/silent
red = 254

allpars = [-1] ;; a 2D array of all parameters for the covariance matrix
allparndx = [-1] ;; an map from the median values of all parameters to the covariance parameters
bestpars = [-1] ;; the best-fit (most likely) values to overplot on the PDF and covariance plots
parnames = [''] ;; parameter names for the axis labels
medianpars = [-1,-1,-1] ;; median parameters and 68% confidence intervals

minchi2 = min(*ss.chi2,bestndx)

if n_elements(csvfile) ne 0 then begin
   openw, csvlun, csvfile, /get_lun
   printf, csvlun, '#parname, median value, upper errorbar, lower errorbar'
endif

!p.multi=[0,2,4] ;; 8 to a page
page = 1
for i=0, n_tags(ss)-1 do begin
   for j=0, n_elements(ss.(i))-1 do begin
      for k=0, n_tags(ss.(i)[j])-1 do begin

         derive = 0B
         fit = 0B
         m=-1

         ;; this captures the detrending variables
         if (size(ss.(i)[j].(k)))[1] eq 10 then begin ;; if it's a pointer
            if ptr_valid(ss.(i)[j].(k)) then begin ;; if it's not empty
               for l=0L, n_tags(*(ss.(i)[j].(k)))-1 do begin ;; loop through each tag
                  if (size((*(ss.(i)[j].(k))).(l)))[2] eq 8 then begin ;; if it's an array of structures
                     for m=0L, n_elements((*(ss.(i)[j].(k))).(l))-1 do begin ;; loop through each structure
                        if tag_exist((*(ss.(i)[j].(k))).(l)[m],'derive') then begin ;; if it's a parameter
                           if (*(ss.(i)[j].(k))).(l)[m].derive or (*(ss.(i)[j].(k))).(l)[m].fit then begin ;; if it's requested to be fit or 
                              
                              ;; remove the burn-in, discard bad chains
                              pars0 = (*(ss.(i)[j].(k))).(l)[m].value
                              if n_elements(mask) ne 0 then pars0[mask] = !values.d_nan
                              pars = (reform(pars0,nsteps,ss.nchains))[burnndx:*,goodchains]

                              ;; GDL doesn't support multi-page PS files
                              if runninggdl and !p.multi[0] eq 0 then begin
                                 device, /close
                                 pdfname0 = file_dirname(pdfname) + path_sep() + file_basename(pdfname,'.ps') + '.' + string(page,format='(i03)') + '.ps'
                                 device, filename=pdfname0, set_font='Times',/tt_font, xsize=xsize,ysize=ysize,/color,bits=24
                                 page += 1
                              endif

                              ;; plot the histogram, get the median and 68% confidence interval
                              summarizepar, pars, label= (*(ss.(i)[j].(k))).(l)[m].label,$
                                            unit=(*(ss.(i)[j].(k))).(l)[m].unit,$
                                            latex=(*(ss.(i)[j].(k))).(l)[m].latex,$
                                            best=(*(ss.(i)[j].(k))).(l)[m].value[bestndx],$
                                            logname=logname,value=value,errlo=errlo,errhi=errhi,$
                                            medianpars=medianpars,charsize=charsize, scinote=scinote
                              
                              ;; populate the median value and error bars
                              (*ss.(i)[j].(k)).(l)[m].medvalue = strtrim(value,2)
                              (*ss.(i)[j].(k)).(l)[m].upper = strtrim(errhi,2)
                              (*ss.(i)[j].(k)).(l)[m].lower = strtrim(errlo,2)
                              (*ss.(i)[j].(k)).(l)[m].scinote = scinote
                              (*ss.(i)[j].(k)).(l)[m].best = (*(ss.(i)[j].(k))).(l)[m].value[bestndx]
                              
                              if n_elements(csvfile) ne 0 then begin
                                 parname = (*(ss.(i)[j].(k))).(l)[m].label
;                                 parname = (*(ss.(i)[j].(k))).(l)[m].latex
                                 parname += '_' + strtrim(j,2) + '_' + strtrim(m,2)
                                 if scinote eq '' then begin
                                    printf, csvlun, parname, strtrim(value,2), strtrim(errhi,2), strtrim(errlo,2),format='(3(a,","),a)'
                                 endif else begin
                                    printf, csvlun, parname, strtrim(value,2), strtrim(errhi,2), strtrim(errlo,2),scinote,format='(4(a,","),a)'
                                 endelse
                              endif

                              ;; store these for the covariance plot
                              if (*(ss.(i)[j].(k))).(l)[m].fit or keyword_set(useallpars) then begin
                                 sz = size(pars)
                                 if n_elements(allpars) eq 1 then allpars = transpose(reform(pars,sz[1]*sz[2])) $
                                 else allpars = [allpars,transpose(reform(pars,sz[1]*sz[2]))]
                                 parnames = [parnames,(*(ss.(i)[j].(k))).(l)[m].latex]

                                 if allparndx[0] eq -1 then allparndx = [n_elements(medianpars[0,*])-2] $
                                 else allparndx = [allparndx,n_elements(medianpars[0,*])-2]
                                 bestpars = [bestpars,(*(ss.(i)[j].(k))).(l)[m].value[bestndx]]
                              endif
                     
                           endif
                        endif
                     endfor
                  endif
               endfor
            endif            
         endif else if n_tags(ss.(i)[j].(k)) ne 0 then begin
            ;; and all other parameters
            if n_tags(ss.(i)[j].(k)) ne 0 then begin
               if tag_exist(ss.(i)[j].(k),'derive') then begin
                  if ss.(i)[j].(k).derive or ss.(i)[j].(k).fit then begin

                     ;; remove the burn-in, discard bad chains
                     pars0 = ss.(i)[j].(k).value
                     if n_elements(mask) ne 0 then pars0[mask] = !values.d_nan
                     pars = (reform(pars0,nsteps,ss.nchains))[burnndx:*,goodchains]

                     ;; GDL doesn't support multi-page PS files
                     if runninggdl and !p.multi[0] eq 0 then begin
                        device, /close
                        pdfname0 = file_dirname(pdfname) + path_sep() + file_basename(pdfname,'.ps') + '.' + string(page,format='(i03)') + '.ps'
                        device, filename=pdfname0, set_font='Times',/tt_font, xsize=xsize,ysize=ysize,/color,bits=24
                        page += 1
                     endif

                     ;; plot the histogram, get the median and 68% confidence interval
                     summarizepar, pars, label=ss.(i)[j].(k).label,$ 
                                   unit=ss.(i)[j].(k).unit,$
                                   latex=ss.(i)[j].(k).latex,$
                                   best=ss.(i)[j].(k).value[bestndx],$
                                   logname=logname,value=value,errlo=errlo,errhi=errhi,$
                                   medianpars=medianpars,charsize=charsize, scinote=scinote
                                   
                     ;; populate the median value and error bars
                     ss.(i)[j].(k).medvalue = strtrim(value,2)
                     ss.(i)[j].(k).upper = strtrim(errhi,2)
                     ss.(i)[j].(k).lower = strtrim(errlo,2)
                     ss.(i)[j].(k).scinote = scinote
                     ss.(i)[j].(k).best = ss.(i)[j].(k).value[bestndx]

                     if n_elements(csvfile) ne 0 then begin
                        parname = ss.(i)[j].(k).label
;                        parname = ss.(i)[j].(k).latex
                        parname += '_' + strtrim(j,2)
                        if scinote eq '' then begin
                           printf, csvlun, parname, strtrim(value,2), strtrim(errhi,2), strtrim(errlo,2),format='(3(a,","),a)'
                        endif else begin
                           printf, csvlun, parname, strtrim(value,2), strtrim(errhi,2), strtrim(errlo,2),scinote,format='(4(a,","),a)'
                        endelse
                     endif

                     ;; store these for the covariance plot
                     if ss.(i)[j].(k).fit or keyword_set(useallpars) then begin
                        sz = size(pars)
                        if n_elements(allpars) eq 1 then allpars = transpose(reform(pars,sz[1]*sz[2])) $
                        else allpars = [allpars,transpose(reform(pars,sz[1]*sz[2]))]
                        parnames = [parnames,ss.(i)[j].(k).latex]
                        if allparndx[0] eq -1 then allparndx = [n_elements(medianpars[0,*])-2] $
                        else allparndx = [allparndx,n_elements(medianpars[0,*])-2]
                        bestpars = [bestpars,ss.(i)[j].(k).value[bestndx]]
                     endif

                  endif
               endif
            endif
         endif
      endfor
   endfor
endfor
device, /close

if n_elements(csvfile) ne 0 then free_lun, csvlun

bestpars = bestpars[1:n_elements(bestpars)-1] ;; the best-fit (most likely) values to overplot on the PDF and covariance plots
parnames = parnames[1:n_elements(parnames)-1] ;; parameter names for the axis labels
medianpars = medianpars[*,1:n_elements(medianpars[0,*])-1] ;; median parameters and 68% confidence intervals

if runninggdl then nocovar = 1

;; if covariance plots aren't wanted (they take a while), we're done
if keyword_set(nocovar) then begin
   !p.font=-1
   !p.multi=0
   set_plot, mydevice
   return
endif

device, filename = covarname
;corr = correlate(allpars)
npars = n_elements(allpars[*,0])

!p.multi=0
charsize = 4d0/npars
exofast_multiplot, [npars,npars], /square, mXtitle='', mYtitle='', /rowmajor
for i=0, npars-1 do begin
  
   ;; x range of the plots/histograms
   xmax = (medianpars[0,allparndx[i]] + 4*medianpars[1,allparndx[i]]) < max(allpars[i,*],/nan)
   xmin = (medianpars[0,allparndx[i]] - 4*medianpars[2,allparndx[i]]) > min(allpars[i,*],/nan)

   for j=0,npars-1 do begin

      ;; the axis labels
      if i eq 0 then ytitle='!3' + exofast_textoidl(parnames[j]) else ytitle = ''
      if j eq npars-1 then xtitle='!3' + exofast_textoidl(parnames[i]) else xtitle = ''

      ;; skip the other half of the corner
      if j lt i+1 then begin

         ;; make a histogram in the diagonal elements
         if j eq i then begin
            hist = histogram(allpars[i,*],locations=x,nbin=100,min=xmin,max=xmax)
            plot, x, hist/total(hist), psym=10,xminor=1, yminor=1,$ 
                  xtickname=replicate(" ",60),ytickname=replicate(" ",60),xtitle=xtitle,$
                  charsize=charsize,xrange=[xmin,xmax],/xstyle
            
            ;; plot lines denoting the median and 68% confidence interval
            median = medianpars[0,allparndx[i]]
            losig = medianpars[0,allparndx[i]] - medianpars[1,allparndx[i]]
            hisig = medianpars[0,allparndx[i]] + medianpars[1,allparndx[i]]
            oplot, [median,median],[-9d9,9d9]
            oplot, [losig,losig],[-9d9,9d9],linestyle=1
            oplot, [hisig,hisig],[-9d9,9d9],linestyle=1

         endif

         exofast_multiplot ;; next plot
         continue
      endif
      
      ;; y range of the plots/histograms
      ymax = (medianpars[0,allparndx[j]] + 4*medianpars[1,allparndx[j]]) < max(allpars[j,*],/nan)
      ymin = (medianpars[0,allparndx[j]] - 4*medianpars[2,allparndx[j]]) > min(allpars[j,*],/nan)
      
      if xmin ne xmax and ymin ne ymax then begin

         xpath = -1
         ypath = -1

         ;; get the paths of the 68% and 95% contours
         exofast_errell, transpose(allpars[i,*]),transpose(allpars[j,*]),xpath=xpath,ypath=ypath,$
                         prob=probs,nxbin=100, nybin=100,$
                         xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,logname=logname
         
         ;; plot the contours if it worked
         if xpath[0] ne -1 and ypath[0] ne -1 then begin

            plot, xpath, ypath, xtitle=xtitle,ytitle=ytitle,charsize=charsize, $
                  xrange=[xmin,xmax],yrange=[ymin,ymax],/ystyle,/xstyle,xminor=1, yminor=1,$ 
                  xtickname=replicate(" ",60),ytickname=replicate(" ",60)
         
            ;; overplot the best fit point
            if n_elements(bestpars) eq npars then begin
               plotsym, 0, 3.5d0/npars, /fill
               oplot, [medianpars[0,allparndx[i]]],[medianpars[0,allparndx[j]]], psym=8
            endif

         endif
      endif
      exofast_multiplot ;; next plot
   endfor

endfor

;; reset back to default plotting parameters
exofast_multiplot, /reset
exofast_multiplot, /default 

;; clean up the postscript device
device, /close
!p.font=-1
!p.multi=0
set_plot, mydevice

end
