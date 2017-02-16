pro 3D MAP

;IMPORT THE .TIF FILES IN IDL ENVIRONMENT
;
   
   ;COHERENCE
path1='/media/marco/marco1/NAPA/Cosismico/Napa_IDL_code/Tiff_pre/Napa_Coherence.tif'
img1_pos=[0]
ENVI_OPEN_FILE, path1, R_FID=fid_img_coer, /NO_REALIZE
ENVI_FILE_QUERY, fid_img_coer, dims=dims_img, ns=ns, nl=nl, xstart=xs, ystart=ys
mapinfo_1=ENVI_GET_MAP_INFO(fid=fid_img_coer, UNDEFINED=mapok)
Coherence=envi_get_data(fid=fid_img_coer, dims=dims_img, pos=[img1_pos])
ENVI_WRITE_ENVI_FILE, Coherence, bnames=['coerenza'],map_info=infomap_img1, 	/IN_MEMORY;
;SAVE, data_coerenza, FILENAME = 'C:\AAA_INGV_CNT\data\xxx.sav'
ENVI_WRITE_ENVI_FILE, Coherence, bnames=['coerenza'],map_info=mapinfo_1, OUT_NAME='/media/marco/marco1/NAPA/Cosismico/Napa_IDL_code/CC.tif'
   
   ;GPS_NS path2='/media/marco/marco1/NAPA/Cosismico/Napa_IDL_code/Tiff_pre/Napa_GPS_NS.tif'
img1_pos=[0]
ENVI_OPEN_FILE, path2, R_FID=fid_img_NS_gps, /NO_REALIZE
ENVI_FILE_QUERY, fid_img_NS_gps, dims=dims_img, ns=ns, nl=nl, xstart=xs, ystart=ys
mapinfo_2=ENVI_GET_MAP_INFO(fid=fid_img_NS_gps, UNDEFINED=mapok)
GPS_NS=envi_get_data(fid=fid_img_NS_gps, dims=dims_img, pos=[img1_pos])
ENVI_WRITE_ENVI_FILE, GPS_NS, bnames=['NS_gps'],map_info=infomap_img1, /IN_MEMORY;
   
   path8='/media/marco/marco1/NAPA/Cosismico/Napa_IDL_code/Tiff_pre/Napa_GPS_Sigma_NS.tif'
img1_pos=[0]
ENVI_OPEN_FILE, path8, R_FID=fid_img_Sigma_NS_gps, /NO_REALIZE
ENVI_FILE_QUERY, fid_img_Sigma_NS_gps, dims=dims_img, ns=ns, nl=nl, xstart=xs, 	ystart=ys
mapinfo_8=ENVI_GET_MAP_INFO(fid=fid_img_Sigma_NS_gps, UNDEFINED=mapok)
GPS_Sigma_NS=envi_get_data(fid=fid_img_Sigma_NS_gps, dims=dims_img, pos=[img1_pos])
ENVI_WRITE_ENVI_FILE, GPS_Sigma_NS,bnames=['Sigma_NS_gps'],map_info=infomap_img1, 	/IN_MEMORY;
   
   ;GPS_EW
path3='/media/marco/marco1/NAPA/Cosismico/Napa_IDL_code/Tiff_pre/Napa_GPS_EW.tif'
img1_pos=[0]
ENVI_OPEN_FILE, path3, R_FID=fid_img_ew_gps, /NO_REALIZE
ENVI_FILE_QUERY, fid_img_ew_gps, dims=dims_img, ns=ns, nl=nl, xstart=xs, ystart=ys
mapinfo_3=ENVI_GET_MAP_INFO(fid=fid_img_ew_gps, UNDEFINED=mapok)
GPS_EW=envi_get_data(fid=fid_img_ew_gps, dims=dims_img, pos=[img1_pos])
ENVI_WRITE_ENVI_FILE, GPS_EW, bnames=['NS_gps'],map_info=infomap_img1, /IN_MEMORY;
   
   path9='/media/marco/marco1/NAPA/Cosismico/Napa_IDL_code/Tiff_pre/Napa_GPS_Sigma_EW.tif'
img1_pos=[0]
ENVI_OPEN_FILE, path9, R_FID=fid_img_Sigma_EW_gps, /NO_REALIZE
ENVI_FILE_QUERY, fid_img_Sigma_EW_gps, dims=dims_img,ns=ns,nl=nl,xstart=xs,ystart=ys
mapinfo_9=ENVI_GET_MAP_INFO(fid=fid_img_Sigma_EW_gps, UNDEFINED=mapok)
GPS_Sigma_EW=envi_get_data(fid=fid_img_Sigma_EW_gps, dims=dims_img, pos=[img1_pos])
ENVI_WRITE_ENVI_FILE, GPS_Sigma_EW,bnames=['Sigma_EW_gps'],map_info=infomap_img1, 	/IN_MEMORY;
   
   ;GPS_UP
path4='/media/marco/marco1/NAPA/Cosismico/Napa_IDL_code/Tiff_pre/Napa_GPS_UP.tif'
img1_pos=[0]
ENVI_OPEN_FILE, path4, R_FID=fid_img_UP_gps, /NO_REALIZE
ENVI_FILE_QUERY, fid_img_UP_gps, dims=dims_img, ns=ns, nl=nl, xstart=xs, ystart=ys
mapinfo_4=ENVI_GET_MAP_INFO(fid=fid_img_UP_gps, UNDEFINED=mapok)
GPS_UP=envi_get_data(fid=fid_img_UP_gps, dims=dims_img, pos=[img1_pos])
ENVI_WRITE_ENVI_FILE, GPS_UP, bnames=['UP_gps'],map_info=infomap_img1, /IN_MEMORY;
   
   path10='/media/marco/marco1/NAPA/Cosismico/Napa_IDL_code/Tiff_pre/Napa_GPS_Sigma_UP.tif'
img1_pos=[0]
ENVI_OPEN_FILE, path10, R_FID=fid_img_Sigma_UP_gps, /NO_REALIZE
ENVI_FILE_QUERY, fid_img_Sigma_UP_gps, dims=dims_img, ns=ns, nl=nl, xstart=xs, 	ystart=ys
mapinfo_10=ENVI_GET_MAP_INFO(fid=fid_img_Sigma_UP_gps, UNDEFINED=mapok)   GPS_Sigma_UP=envi_get_data(fid=fid_img_Sigma_UP_gps, dims=dims_img, pos=[img1_pos])
ENVI_WRITE_ENVI_FILE, GPS_Sigma_UP,bnames=['Sigma_UP_gps'],map_info=infomap_img1, 	/IN_MEMORY;
   
   ;InSAR
path5='/media/marco/marco1/NAPA/Cosismico/Napa_IDL_code/Tiff_pre/Napa_InSAR.tif'
img1_pos=[0]
ENVI_OPEN_FILE, path5, R_FID=fid_img_insar, /NO_REALIZE
ENVI_FILE_QUERY, fid_img_insar, dims=dims_img, ns=ns, nl=nl, xstart=xs, ystart=ys
mapinfo_5=ENVI_GET_MAP_INFO(fid=fid_img_insar, UNDEFINED=mapok)
D_InSAR_LOS=envi_get_data(fid=fid_img_insar, dims=dims_img, pos=[img1_pos])
ENVI_WRITE_ENVI_FILE, D_InSAR_LOS, bnames=['InSAR'],map_info=infomap_img1, 	/IN_MEMORY;
;SAVE, data_coerenza, FILENAME = 'C:\AAA_INGV_CNT\data\xxx.sav'
   
   ;MAI
path6='/media/marco/marco1/NAPA/Cosismico/Napa_IDL_code/Tiff_pre/Napa_MAI.tif'
img1_pos=[0]
ENVI_OPEN_FILE, path6, R_FID=fid_img_mai1, /NO_REALIZE
ENVI_FILE_QUERY, fid_img_mai1, dims=dims_img, ns=ns, nl=nl, xstart=xs, ystart=ys
mapinfo_6=ENVI_GET_MAP_INFO(fid=fid_img_mai1, UNDEFINED=mapok)
MAI_AZ=envi_get_data(fid=fid_img_mai1, dims=dims_img, pos=[img1_pos])
ENVI_WRITE_ENVI_FILE, MAI_AZ, bnames=['MAI_AZ'],map_info=infomap_img1, /IN_MEMORY;
;SAVE, data_coerenza, FILENAME = 'C:\AAA_INGV_CNT\data\xxx.sav'
   
   ;POT
path7='/media/marco/marco1/NAPA/Cosismico/Napa_IDL_code/Tiff_pre/NAPA_POT.tif'
img1_pos=[0]
ENVI_OPEN_FILE, path7, R_FID=fid_img_pot1, /NO_REALIZE
ENVI_FILE_QUERY, fid_img_pot1, dims=dims_img, ns=ns, nl=nl, xstart=xs, ystart=ys
mapinfo_7=ENVI_GET_MAP_INFO(fid=fid_img_pot1, UNDEFINED=mapok)
POT_AZ=envi_get_data(fid=fid_img_pot1, dims=dims_img, pos=[img1_pos])
ENVI_WRITE_ENVI_FILE, POT_AZ, bnames=['POT_AZ'],map_info=infomap_img1, /IN_MEMORY;
;SAVE, data_coerenza, FILENAME = 'C:\AAA_INGV_CNT\data\xxx.sav'
   
   
   ;REDEFINING NAN
for i=0, 1815 do begin
  for j=0, 1168 do begin
    if (D_InSAR_LOS(i,j) EQ 0.) then begin
      print, 'Reina'
      D_InSAR_LOS(i,j)=!VALUES.F_NAN
    endif  else begin
      D_InSAR_LOS(i,j)=D_InSAR_LOS(i,j)
      print, 'Hysaj'
    endelse
  endfor
endfor

   
for i=0, 1815 do begin
  for j=0, 1168 do begin
    if (MAI_AZ(i,j) EQ 0.) then begin
      print, 'Albiol'
      MAI_AZ(i,j)=!VALUES.F_NAN
    endif else begin
      MAI_AZ(i,j)=MAI_AZ(i,j)
      print, 'Koulibaly'
    endelse
  endfor
endfor

   
for i=0, 1815 do begin
  for j=0, 1168 do begin
    if (POT_AZ(i,j) EQ 0.) then begin
      print, 'Ghoulam'
      POT_AZ(i,j)=!VALUES.F_NAN
    endif else begin
      POT_AZ(i,j)=POT_AZ(i,j)
      print, 'Allan'
    endelse
  endfor
endfor

   
for i=0, 1815 do begin
  for j=0, 1168 do begin
    if (GPS_NS(i,j) LT -1e+15) then begin
      print, 'Jorginho'
      GPS_NS(i,j)=!values.f_nan
    endif  else begin
      GPS_NS(i,j)=GPS_NS(i,j)
      print, 'Hamsik'
    endelse
  endfor
endfor

  
for i=0, 1815 do begin
  for j=0, 1168 do begin
    if (GPS_UP(i,j) LT -1e+15) then begin
      print, 'Callejon'
      GPS_UP(i,j)=!values.f_nan
    endif  else begin
      GPS_UP(i,j)=GPS_UP(i,j)
      print, 'Higuain'
    endelse
  endfor
endfor

   
for i=0, 1815 do begin
  for j=0, 1168 do begin
    if (GPS_EW(i,j) LT -1e+15) then begin
      print, 'Insigne'
      GPS_EW(i,j)=!values.f_nan
    endif  else begin
      GPS_EW(i,j)=GPS_EW(i,j)
      print, 'Gabbiadini'
    endelse
  endfor
endfor

   
for i=0, 1815 do begin
  for j=0, 1168 do begin
    if (GPS_Sigma_NS(i,j) LT -1e+15) then begin
      print, 'Mertens'
      GPS_Sigma_NS(i,j)=!values.f_nan
    endif  else begin
      GPS_Sigma_NS(i,j)=GPS_Sigma_NS(i,j)
      print, 'Lopez'
    endelse
  endfor
endfor

   
for i=0, 1815 do begin
  for j=0, 1168 do begin
    if (GPS_Sigma_EW(i,j) LT -1e+15) then begin
      print, 'Chiriches'
      GPS_Sigma_EW(i,j)=!values.f_nan
    endif  else begin
      GPS_Sigma_EW(i,j)=GPS_Sigma_EW(i,j)
      print, 'Valdifiori'
    endelse
  endfor
endfor

   
for i=0, 1815 do begin
  for j=0, 1168 do begin
    if (GPS_Sigma_UP(i,j) LT -1e+15) then begin
    print, 'Strinic'
      GPS_Sigma_UP(i,j)=!values.f_nan
    endif  else begin
      GPS_Sigma_UP(i,j)=GPS_Sigma_UP(i,j)
      print, 'Maggio'
    endelse
  endfor
endfor
  
   

;SET MAI AND POT OUTPUTS AS -1*MEASURED DISPLACEMENT SINCE ALONG DESCENDING ORBIT THE MOVEMENT IS CONSIDERED POSITIVE SOUTHWARD (I.E. IN THE SAME DIRECTION OF THE SATELLITE FLIGHT)
D_MAI_AZ=-1*MAI_AZ
D_POT_AZ=-1*POT_AZ

;SET SIGMA_GPS AS SIGMA_INTERPOLATION + MEAN SIGMA MEASURED AND THEN CONVERT IN METERS
Sigma_GPS_NS=((GPS_Sigma_NS+1.3351807229)/1000)
Sigma_GPS_EW=((GPS_Sigma_EW+1.180753012)/1000)
Sigma_GPS_UP=((GPS_Sigma_UP+4.2909717868)/1000)

;CONVERT IN METERS THE DISPLACEMENT MEASURED BY GPS
D_GPS_NS=GPS_NS/1000
D_GPS_EW=GPS_EW/1000
D_GPS_UP=GPS_UP/1000

;SATELLITE AND PROCESSING PARAMETERS
cos_alpha=-0.97334353835
sen_alpha=-0.22935203583
cos_teta=0.91922057496
sen_teta=0.39374298033

NL=225 ;Numero di Looks
n=0.5 ;Normalized Squint
l=12.3 ;Antenna azimuth lenght
lambda=0.05550415767 ;wavelenght

;INITIALIZE THE MATRICES
Sigma_MAI=fltarr(1816,1169)
Sigma_POT=fltarr(1816,1169)
Sigma_InSAR=fltarr(1816,1169)

d_EW=fltarr(1816,1169)
d_NS=fltarr(1816,1169)
d_UP=fltarr(1816,1169)

A=fltarr(1816,1169)
B=fltarr(1816,1169)
C=fltarr(1816,1169)
D=fltarr(1816,1169)
E=fltarr(1816,1169)
F=fltarr(1816,1169)
G=fltarr(1816,1169)
H=fltarr(1816,1169)
P=fltarr(1816,1169)
M=fltarr(1816,1169)
;prova=make_array(1816,1169, value=!Values.F_NAN )


;ESTIMATION OF THE 3D OPTIMIZED DISPLACEMENT COMPONENTS


;SIGMA
for i=0, 1815 do begin
  for j=0, 1168 do begin
    Sigma_MAI(i,j)=((l/(4.*!pi*n))*(1./(sqrt(NL)))*(sqrt(1.-	Coherence(i,j))^2)/Coherence(i,j)))
    Sigma_POT(i,j)=((l/2.)*(sqrt(3./(10.*NL)))*(sqrt(2.+(5.*(Coherence(i,j)^2))-	(7.*(Coherence(i,j)^4)))/!pi*(Coherence(i,j)^2)))
    Sigma_InSAR(i,j)=(((lambda/(4.*!pi))*(1./sqrt(2.*NL))*(sqrt(1.-	(Coherence(i,j))^2)/Coherence(i,j))))

;SUPPORT VARIABLES
    A(i,j)=(((-1./Sigma_InSAR(i,j)^2)*cos_alpha^2*cos_teta^2)+((1./	(Sigma_MAI(i,j))^2)*sen_alpha^2)+		((1./Sigma_POT(i,j)^2)*sen_alpha^2)+(1./Sigma_GPS_EW(i,j)^2))
    B(i,j)=(1./Sigma_InSAR(i,j)^2*cos_alpha*cos_teta*sen_teta)
    C(i,j)=((1./Sigma_MAI(i,j)^2)*cos_alpha*sen_alpha+	((1./Sigma_POT(i,j)^2)*cos_alpha*sen_alpha))
    D(i,j)=(((1./Sigma_InSAR(i,j)^2)*D_InSAR_LOS(i,j)*cos_alpha*sen_teta)+	((1./Sigma_MAI(i,j)^2)*MAI_AZ(i,j)*sen_alpha)+	((1./Sigma_POT(i,j)^2)*POT_AZ(i,j)*sen_alpha)+	((1./Sigma_GPS_EW(i,j)^2)*D_GPS_EW(i,j)))
    E(i,j)=(((1./Sigma_MAI(i,j)^2)*cos_alpha^2)+	((1./Sigma_POT(i,j)^2)*cos_alpha^2)+(1./Sigma_GPS_NS(i,j)^2))
    F(i,j)=(((1./Sigma_MAI(i,j)^2)*sen_alpha*cos_alpha)-	((1./Sigma_POT(i,j)^2)*sen_alpha*cos_alpha))
    G(i,j)=(((1./Sigma_MAI(i,j)^2)*MAI_AZ(i,j)*cos_alpha)+	((1./Sigma_POT(i,j)^2)*POT_AZ(i,j)*cos_alpha)+	((1./Sigma_GPS_NS(i,j)^2)*D_GPS_NS(i,j)))
    H(i,j)=(((1./Sigma_InSAR(i,j)^2)*cos_teta^2)+(1./Sigma_GPS_UP(i,j)^2))
    P(i,j)=((1./Sigma_InSAR(i,j)^2)*cos_alpha*sen_teta*cos_teta)
    M(i,j)=(((1./Sigma_InSAR(i,j)^2)*D_InSAR_LOS(i,j)*cos_teta)+	((1./Sigma_GPS_UP(i,j)^2)*D_GPS_UP(i,j)))

;3D COMPONENTS    
    d_EW(i,j)=(((-B(i,j)/A(i,j))*(M(i,j)/(H(i,j)))-((C(i,j)/A(i,j))*(G(i,j)/E(i,j)))	+(D(i,j)/A(i,j)))/(1+((B(i,j)/A(i,j))*(P(i,j)/H(i,j)))+	((C(i,j)/A(i,j))*(F(i,j)/E(i,j)))))
    d_NS(i,j)=((F(i,j)/E(i,j))*d_EW(i,j))+(G(i,j)/E(i,j))
    d_UP(i,j)=(((P(i,j))/(H(i,j))*d_EW(i,j))+(M(i,j)/H(i,j)))
  
  endfor
endfor

;RIDEFINING NAN IN SIGMA MAPS
for i=0, 1815 do begin
  for j=0, 1168 do begin
    if (Sigma_MAI(i,j) LT -1e+15) then begin
      print, 'ciao'
      Sigma_MAI(i,j)=!values.f_nan
    endif  else begin
      Sigma_MAI(i,j)=Sigma_MAI(i,j)
      print, 'pluto'
    endelse
  endfor
endfor

for i=0, 1815 do begin
  for j=0, 1168 do begin
    if (Sigma_MAI(i,j) GT 1e+15) then begin
      print, 'ciao'
      Sigma_MAI(i,j)=!values.f_nan
    endif  else begin
      Sigma_MAI(i,j)=Sigma_MAI(i,j)
      print, 'pluto'
    endelse
  endfor
endfor

for i=0, 1815 do begin
  for j=0, 1168 do begin
    if (Sigma_InSAR(i,j) LT -1e+15) then begin
      print, 'ciao'
      Sigma_InSAR(i,j)=!values.f_nan
    endif  else begin
      Sigma_InSAR(i,j)=Sigma_InSAR(i,j)
      print, 'pluto'
    endelse
  endfor
endfor

for i=0, 1815 do begin
  for j=0, 1168 do begin
    if (Sigma_InSAR(i,j) GT 1e+15) then begin
      print, 'ciao'
      Sigma_InSAR(i,j)=!values.f_nan
    endif  else begin
      Sigma_InSAR(i,j)=Sigma_InSAR(i,j)
      print, 'pluto'
    endelse
  endfor
endfor

for i=0, 1815 do begin
  for j=0, 1168 do begin
    if (Sigma_POT(i,j) LT -1e+15) then begin
      print, 'ciao'
      Sigma_POT(i,j)=!values.f_nan
    endif  else begin
      Sigma_POT(i,j)=Sigma_POT(i,j)
      print, 'pluto'
    endelse
  endfor
endfor

for i=0, 1815 do begin
  for j=0, 1168 do begin
    if (Sigma_POT(i,j) GT 1e+15) then begin
      print, 'ciao'
      Sigma_POT(i,j)=!values.f_nan
    endif  else begin
      Sigma_POT(i,j)=Sigma_POT(i,j)
      print, 'pluto'
    endelse
  endfor
endfor

;H=(Sigma_POT+Sigma_GPS_NS)/2
;F=(Sigma_MAI+Sigma_GPS_EW+Sigma_InSAR+Sigma_POT)/3.5

;EXPORT TO .TIF FORMAT
ENVI_WRITE_ENVI_FILE, Sigma_GPS_EW, BNAMES=['PROVA'],MAP_INFO=MAPINFO_1, OUT_NAME='/media/marco/marco1/NAPA/Cosismico/Napa_IDL_code/Output/SIGMA_GPS_EW.tif'

ENVI_WRITE_ENVI_FILE, Sigma_GPS_UP, BNAMES=['PROVA'],MAP_INFO=MAPINFO_1, OUT_NAME='/media/marco/marco1/NAPA/Cosismico/Napa_IDL_code/Output/SIGMA_GPS_UP.tif'

ENVI_WRITE_ENVI_FILE, Sigma_GPS_NS, bnames=['prova'],map_info=mapinfo_1, OUT_NAME='/media/marco/marco1/NAPA/Cosismico/Napa_IDL_code/Output/SIGMA_GPS_NS.tif'

ENVI_WRITE_ENVI_FILE, D_MAI_AZ, bnames=['prova'],map_info=mapinfo_1, OUT_NAME='/media/marco/marco1/NAPA/Cosismico/Napa_IDL_code/Output/MAI_AZ.tif'

ENVI_WRITE_ENVI_FILE, D_InSAR_LOS, bnames=['prova'],map_info=mapinfo_1, OUT_NAME='/media/marco/marco1/NAPA/Cosismico/Napa_IDL_code/Output/InSAR_LOS.tif'

ENVI_WRITE_ENVI_FILE, D_POT_AZ, bnames=['prova'],map_info=mapinfo_1, OUT_NAME='/media/marco/marco1/NAPA/Cosismico/Napa_IDL_code/Output/POT_AZ.tif'

ENVI_WRITE_ENVI_FILE, Sigma_MAI, bnames=['prova'],map_info=mapinfo_1, OUT_NAME='/media/marco/marco1/NAPA/Cosismico/Napa_IDL_code/Output/SIGMA_MAI.tif'

ENVI_WRITE_ENVI_FILE, Sigma_InSAR, bnames=['prova'],map_info=mapinfo_1, OUT_NAME='/media/marco/marco1/NAPA/Cosismico/Napa_IDL_code/Output/SIGMA_InSAR.tif'

ENVI_WRITE_ENVI_FILE, Sigma_POT, bnames=['prova'],map_info=mapinfo_1, OUT_NAME='/media/marco/marco1/NAPA/Cosismico/Napa_IDL_code/Output/SIGMA_POT.tif'

ENVI_WRITE_ENVI_FILE, H, bnames=['prova'],map_info=mapinfo_1, OUT_NAME='/media/marco/marco1/NAPA/Cosismico/Napa_IDL_code/Output/SIGMA_NS.tif'

ENVI_WRITE_ENVI_FILE, F, bnames=['prova'],map_info=mapinfo_1, OUT_NAME='/media/marco/marco1/NAPA/Cosismico/Napa_IDL_code/Output/SIGMA_EW.tif'

;EXPORT ALSO THE 3D OPTIMIZED COMPONENTS IN .TIF FORMAT
ENVI_WRITE_ENVI_FILE, d_EW, bnames=['d_EW'],map_info=mapinfo_1, OUT_NAME='/media/marco/marco1/NAPA/Cosismico/Napa_IDL_code/Output/d_EW.tif'

ENVI_WRITE_ENVI_FILE, d_NS, bnames=['d_NS'],map_info=mapinfo_1, OUT_NAME='/media/marco/marco1/NAPA/Cosismico/Napa_IDL_code/Output/d_NS.tif'

ENVI_WRITE_ENVI_FILE, d_UP, bnames=['d_UP'],map_info=mapinfo_1, OUT_NAME='/media/marco/marco1/NAPA/Cosismico/Napa_IDL_code/Output/d_UP.tif'

stop

end

