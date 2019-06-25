# glacier_runoff_subset()
# Function to subset PyGEM glacier runoff output to given years
# Project: NASA HiMAT
# Danielle S Grogan
# Last updated 2019-06-25

library(raster)

############################################################################################################
glacier_runoff_subset = function(gl.path,     # path to glacier model output
                                 model,       # If rcp != historical, also supply a GCM model name
                                 rcp,         # rcp = one of: "historical", "rcp45", "rcp85"
                                 st.yr,       # start year to subset
                                 end.yr,      # end year to subset
                                 out.yr = 0){ # 0 for output in m3/month, 1 for output in m3/year

  # glacier runoff from glacier model
  if(rcp == 'historical'){
    glacier.runoff = brick(file.path(gl.path, "ERA-Interim_c2_ba0_200sets_2000_2017_stats_m.nc"),
                           varname='runoff')
  }else{
    glacier.runoff = brick(paste(gl.path, "/", model, "_", rcp, "_c2_ba2_100sets_2000_2100_m.nc", sep=""),
                           varname='runoff')
  }
  
  # subset to defined start and end years
  n = 12*(end.yr - st.yr)
  layer.1 = min(which(grepl(st.yr, c(names(glacier.runoff)))))
  glacier.runoff.sub = subset(glacier.runoff, layer.1:(layer.1 + n-1))
  
  if(out.yr == 0){
    out = glacier.runoff.sub
  }else if(out.yr == 1){
    # sum the months in each year
    ids = unlist(lapply(seq(1,n/12), FUN = function(x) rep(x,12)))
    glacier.runoff.m3year = stackApply(glacier.runoff.sub, indices = ids, fun = sum)
    out = glacier.runoff.m3year
  }
  out
}
############################################################################################################

