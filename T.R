for(f in list.files("R/",full.names =T)){
  print(tools:::showNonASCIIfile(f))
}