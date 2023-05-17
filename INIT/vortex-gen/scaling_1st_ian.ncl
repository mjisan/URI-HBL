;*********************************************Read Input Files*****************************************************
  file_F = "track_file_florence_new"
  file_G = "track_file_florence_new_extended"
  
  fname      = "track_file_florence"
  data       = asciiread(fname,-1,"string")
  col6       = str_get_cols(data,0,60)
  track_file = asciiread("track_file_florence",(/2,19/),"float")
  mws        = track_file(:,9)
  rad        = track_file(:,10)
  lonn       = track_file(:,3)
  latt       = track_file(:,2)

  lonn1      = -lonn/10
  latt1      =  latt/10

  lonn2      = lonn1(0)
  lonn3      = lonn1(1)
  latt2      = latt1(0)
  latt3      = latt1(1)

  f1         = addfile("boundary_model.nc","r")

;*********************************************Read Input Variables*************************************************

  u      = f1->ubot
  v      = f1->vbot
  u1     = f1->utop
  v1     = f1->vtop
  X      = f1->x
  Y      = f1->y

  rad18a = track_file(:,11)
  rad18b = track_file(:,12)
  rad18c = track_file(:,13)
  rad18d = track_file(:,14)
  
  rad26a = track_file(:,15)
  rad26b = track_file(:,16)
  rad26c = track_file(:,17)
  rad26d = track_file(:,18)

;*********************************************Average Quadrant Wind*************************************************

  r26     = rad26a+rad26b+rad26c+rad26d
  r26_avg = r26/4
  r18     = rad18a+rad18b+rad18c+rad18d
  r18_avg = r18/4

;*********************************************Creating Grid*************************************************

  nlon           = dimsizes(X)
  nlat           = dimsizes(Y)
  latgrid_center =  fspan( (-floor(nlat/2)+1),floor(nlat/2), 1802)
  lat            = latgrid_center*1e3*1e-3
  longrid_center = fspan( (-floor(nlon/2)+1), floor(nlon/2), 1802)
  lon            = longrid_center*1e3*1e-3

;*******************************************Wind Speed and Translation Speed Calculator*****************************

  ws_bot  = sqrt(u^2+v^2)
  ws_bot1 = ws_bot(15,:,:)

  ws_top  = sqrt(u1^2+v1^2)
  ws_top1 = ws_top(15,:,:)
  dist    = gc_latlon(latt3,lonn3,latt2,lonn2,2,3)
  wstr1   = dist/(6*3600)


;******************************************Assign Dimensions********************************************************

  lat!0       = "lat"              
  lat@units   = "meters"    
  lat&lat     = lat         
  lon!0       = "lon"
  lon@units   = "meters"
  lon&lon     = lon
  u!0         = "time"
  u!1         = "lat"
  u!2         = "lon"
  u&lat       = lat
  u&lon       = lon
  v!0         = "time"
  v!1         = "lat"
  v!2         = "lon"
  v&lat       = lat
  v&lon       = lon
  ws_bot1!0   = "lat"
  ws_bot1!1   = "lon"
  ws_bot1&lat = lat
  ws_bot1&lon = lon
  ws_top1!0   = "lat"
  ws_top1!1   = "lon"
  ws_top1&lat = lat
  ws_top1&lon = lon

  
;***********************************************Identify grid indices***********************************************

  r181    = floattointeger(r18_avg)
  r261    = floattointeger(r26_avg)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;Calulate Ratio of Wind and Translation;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  vg10    = max(ws_bot1)
  print(vg10)
  
  R_Vg    = mws/vg10
  print(R_Vg)
  
  R_Ut    = wstr1/R_Vg

  print(wstr1)
  print(R_Ut)
  
  New_Vg1 = mws - R_Ut
  print(New_Vg1)
  
  New_Vg  = New_Vg1*R_Vg
  print(New_Vg)

  New_Vg  = decimalPlaces(New_Vg,0,True)
     vg1  = floattointeger(New_Vg)


  ws18_nhc       = ws_top1({0},{r181})
  ws26_nhc       = ws_top1({0},{r261})
   New_r18       = tofloat(ws18_nhc*R_Vg)
  New_r26        = tofloat(ws26_nhc*R_Vg)
  r18_diff       = 0.0
  r26_diff       = 0.0
  New_RMW        = floattointeger(rad)


  
lines = sprinti("%0.2i", vg1)+ " " + sprinti("%0.3i", New_RMW)+ " " + sprinti("%0.4i", r181) + " " +sprinti("%0.4i", r181) + " "+sprinti("%0.4i", r181) + " "+sprinti("%0.4i", r181) + " "+ sprinti("%0.4i", r261) + " " +sprinti("%0.4i", r261) + " " +sprinti("%0.4i", r261) + " " +sprinti("%0.4i", r261) +" "+sprintf("%0.2f", R_Vg)+" "+sprintf("%0.2f", New_r18)+" "+sprintf("%0.2f", New_r26)+ " " +sprintf("%0.2f", r18_diff)+" "+sprintf("%0.2f", r26_diff) 

lines2 = sprinti("%0.2i", vg1)+ " " + sprinti("%0.3i", New_RMW)+ " " + sprinti("%0.4i", r181) + " " +sprinti("%0.4i", r181) + " "+sprinti("%0.4i", r181) + " "+sprinti("%0.4i", r181) + " "+ sprinti("%0.4i", r261) + " " +sprinti("%0.4i", r261) + " " +sprinti("%0.4i", r261) + " " +sprinti("%0.4i", r261) +" "+sprintf("%0.2f", R_Vg)+" "+sprintf("%0.2f", New_r18)+" "+sprintf("%0.2f", New_r26)+ " " +sprintf("%0.2f", r18_diff)+" "+sprintf("%0.2f", r26_diff)+" "+wstr1+" "+R_Ut+" "+New_Vg1+" "+ws18_nhc+" "+ws26_nhc




write_table(file_F,"a",[/col6,lines/], \
              "%s%s")

write_table(file_G,"a",[/col6,lines2/], \
              "%s%s")
