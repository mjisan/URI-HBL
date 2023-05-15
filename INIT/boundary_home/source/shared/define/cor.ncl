  fname   = "track_file_florence-sv"
  nrow    = numAsciiRow(fname)
  ncol    = numAsciiCol(fname)
  a       = asciiread(fname,(/nrow,1/),"float")
  lat1    = a(:,0)
  PI      = 3.141592653589793238
  omega   = 7.2921e-5

  lat     = lat1*(PI/180)
  fcor    = 2*omega*sin(lat) 
  aaa     = sprintf("%0.10f",fcor)
  print(fcor)
asciiwrite("cor.txt",aaa)