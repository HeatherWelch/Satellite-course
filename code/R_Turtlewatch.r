######## NOAA SATELLITE COURSE 2016 ###########

####### USING R FOR SPATIAL ANALYSES - A 'TURTLEWATCH' EXAMPLE  ########

###### SET YOUR WORKING DIRECTORY
setwd('./')

###### INSTALL PACKAGES  #### 

## function to install packages needed
pkgTest <- function(x)
{
  if (!require(x,character.only = TRUE))
  {
    install.packages(x,dep=TRUE)
    if(!require(x,character.only = TRUE)) stop("Package not found")
  }
}

### LOAD LIBRARIES
pkgTest("ncdf4")
pkgTest("RCurl")
pkgTest("raster")
pkgTest("colorRamps")
pkgTest("maps")
pkgTest("mapdata")


###### AUTOMATIC DATA DOWNLOAD FROM ERDDAP

## Navigate to https://coastwatch.pfeg.noaa.gov/erddap/index.html
## Find the NASA Global High Resolution Sea Surface Temperature (GHRSST) Daily 1km dataset
## HINT: You can search for 'jplG1SST'
## CLICK the 'data' link in the GridDAP Data column
## Under Grid Variables, uncheck 'mask' and 'analysis_error', but leave 'SST' checked
## Change File Type to .nc - Download a NetCDF-3 binary file with COARDS/CF/ACDD metadata 
## CLICK 'Just generate the URL' button towards the bottom of the page and copy-paste the URL to your text editor
## It should look something like this:
## 'https://coastwatch.pfeg.noaa.gov/erddap/griddap/jplG1SST.nc?SST[(2016-08-22T00:00:00Z):1:(2016-08-22T00:00:00Z)][(-79.995):1:(79.995)][(-179.995):1:(179.995)]'

##### Now we want to refine the spatial and temporal extent to get the data we want.

## DATES - choose 10 days of summer in 2011, 2012, 2013, 2014 OR 2015

dates<-seq(as.Date('2015/06/01'), as.Date('2015/06/11'), by = 'day',format='%Y-%m-%dd')

## Try with one date to start with:
startdate <- dates[1]
enddate <- dates[2]

## SPATIAL EXTENT 
## Choose an area off the coast of California and define the latitude and longitude of the vertices

min.lat <- 25; max.lat <- 40
min.lon <- -130; max.lon <- -115  ## NB. longitude is -180° to 180° for this product, but others use 0-360°

## BUILD THE URL TO ACCESS THE DATA SUBSET  ## NB. will be a different format for different data types

our.url <- paste('https://coastwatch.pfeg.noaa.gov/erddap/griddap/jplG1SST.nc?SST[(',startdate,'T00:00:00Z):1:(',enddate,'T00:00:00Z)][(',min.lat,'):1:(',max.lat,')][(',min.lon,'):1:(',max.lon,')]',sep='')

## SET UP LOCAL NETCDF FILE AND DOWNLOAD DATA INTO IT
f <- CFILE(paste('jplG1SST_',startdate,'.nc',sep=''),mode="wb")
curlPerform(url=our.url,writedata=f@ref,noprogress=FALSE) ## this may take a couple of minutes - the command line cursor will hang
close(f)  ## important not to forget this step


##### Now run the other dates in a 'while' loop

## Function to wait a few seconds before hitting server
waitfor <- function(x){
    p1 <- proc.time()
    Sys.sleep(x)
    print(proc.time() - p1) # The cpu usage should be negligible
}

### WHILE LOOP
i<-1
waitsecs<-2
while (i < length(dates)){
	startdate<-dates[i]
	enddate<-dates[i+1]
	filenm<-paste('jplG1SST_',startdate,'.nc',sep='')
	url<-paste('https://coastwatch.pfeg.noaa.gov/erddap/griddap/jplG1SST.nc?SST[(',startdate,'T00:00:00Z):1:(',enddate,'T00:00:00Z)][(',min.lat,'):1:(',max.lat,')][(',min.lon,'):1:(',max.lon,')]',sep='')
	print(paste('Starting download for',startdate,'to',enddate))
	f = CFILE(filenm,mode="wb")
	curlPerform(url=url,writedata=f@ref,noprogress=FALSE) 
	close(f)
	print(paste('Download for',startdate,'to',enddate,'finished'))
	### in case our server query didn't work, we can check, wait a couple of seconds and restart
	i<-i+1
	if (is.na(file.info(filenm)$size)) {
		i<-i-1
		waitfor(waitsecs)
		waitsecs<-waitsecs+2
	}
	else if (file.info(filenm)$size < 2000){
		i<-i-1
		waitfor(waitsecs)
		waitsecs<-waitsecs+2
	}
	else waitsecs<-2
	if (waitsecs > 90) waitsecs <- 30
}


##### IMPORT NETCDF AS RASTER IN R

## proj4string for projection system (latitude-longitude)
projstring <- '+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0'

## import NetCDF
nc.data <- nc_open(paste('jplG1SST_',startdate,'.nc',sep=''))
print(nc.data)

## get spatial and temporal bounds
lat <- ncvar_get(nc.data,'latitude'); lon <- ncvar_get(nc.data,'longitude')
tim <- ncvar_get(nc.data,'time'); day <- format(as.POSIXlt(tim,origin='1970-01-01',tz='UTC'),'%Y-%m-%d')
latlons <- expand.grid(lon, lat); colnames(latlons) <- c("lon", "lat")
coordinates(latlons) <- ~lon+lat; proj4string(latlons) <- CRS(projstring)

index1 <- which(abs(as.Date(startdate)-as.Date(day)) == min(abs(as.Date(startdate)-as.Date(day))))[1] 

## get variable of interest
varname <- 'SST'
var <-  ncvar_get(nc.data,varname)
z <- var[,,index1]

## make raster from NetCDF - transpose and rotate variable array and set projection and extent
zRot <-  apply(t(z),2,rev); tR <- raster(zRot); projection(tR) <- CRS(projstring)
extent(tR) <- c(xmn=min(latlons$lon),xmx=max(latlons$lon),ymn=min(latlons$lat),ymx=max(latlons$lat))

## a quick plot to check
image(tR, col=matlab.like2(255))
map('worldHires',add=TRUE, col='grey90', fill=TRUE)


##### DEFINING THE THERMAL NICHE OF LOGGERHEAD TURTLES

## Mask all raster cells outside temperature range

min.temp <- 17.5
max.temp <- 18.5

tR2 <- tR
tR2@data@values <- getValues(tR2)
tR2[!values(tR) >= min.temp] <- NA
tR2[!values(tR) <= max.temp] <- NA

map('worldHires',add=TRUE, col='grey90', fill=TRUE)

#### Make plots to output as a png image, and save the raster objects

## output directory
outdir <- './'

## SST over whole domain
png(paste(outdir,'GHRSST_',startdate,'.png',sep=''),width=960,height=960,units='px',pointsize=20)
par(mar=c(3,3,.5,.5),las=1,font.axis=2)
image(tR, col=matlab.like2(255))
map('worldHires',add=TRUE,col=grey(0.7),fill=TRUE)
box(); dev.off()
writeRaster(tR,filename=paste(outdir,'GHRSST_',startdate,'.grd',sep=''),overwrite=TRUE)

## turtle habitat
png(paste(outdir,'TurtleHabitat_',startdate,'.png',sep=''),width=960,height=960,units='px',pointsize=20)
par(mar=c(3,3,.5,.5),las=1,font.axis=2)
image(tR2, col='firebrick4')
map('worldHires',add=TRUE,col=grey(0.7),fill=TRUE)
box(); dev.off()
writeRaster(tR2,filename=paste(outdir,'TurtleHabitat_',startdate,'.grd',sep=''),overwrite=TRUE)

## overlay turtle habitat on SST image
png(paste(outdir,'SST_TurtleHabitat_',startdate,'.png',sep=''),width=960,height=960,units='px',pointsize=20)
par(mar=c(3,3,.5,.5),las=1,font.axis=2)
image(tR, col=matlab.like2(255))
image(tR2, col='firebrick4',add=TRUE)
map('worldHires',add=TRUE,col=grey(0.7),fill=TRUE)
box(); dev.off()

### We can use a function to plot and save automatically over the whole date range

plot.turtle <- function(predDate,outdir,varname, min.var, max.var){
	print(paste('Plotting',predDate))
	nc.data <- nc_open(paste('jplG1SST_',predDate,'.nc',sep=''))
	lat <- ncvar_get(nc.data,'latitude'); lon <- ncvar_get(nc.data,'longitude')
	tim <- ncvar_get(nc.data,'time'); day <- format(as.POSIXlt(tim,origin='1970-01-01',tz='UTC'),'%Y-%m-%d')
	latlons <- expand.grid(lon, lat); colnames(latlons) <- c("lon", "lat")
	coordinates(latlons) <- ~lon+lat; proj4string(latlons) <- CRS(projstring)
	index1 <- which(abs(as.Date(predDate)-as.Date(day)) == min(abs(as.Date(predDate)-as.Date(day))))[1] 
	var <-  ncvar_get(nc.data,varname)
	z <- var[,,index1]
	zRot <-  apply(t(z),2,rev); tR <- raster(zRot); projection(tR) <- CRS(projstring)
	extent(tR) <- c(xmn=min(latlons$lon),xmx=max(latlons$lon),ymn=min(latlons$lat),ymx=max(latlons$lat))
	tR2 <- tR; tR2@data@values <- getValues(tR2)
	tR2[!values(tR) >= min.var] <- NA; tR2[!values(tR) <= max.var] <- NA
	png(paste(outdir,'GHRSST_TurtleHabitat_',predDate,'.png',sep=''),width=960,height=960,units='px',pointsize=20)
	par(mar=c(3,3,.5,.5),las=1,font.axis=2)
	image(tR, col=matlab.like2(255)); image(tR2, col='firebrick4',add=TRUE)
	map('worldHires',add=TRUE,col=grey(0.7),fill=TRUE)
	box(); dev.off()
	writeRaster(tR,filename=paste(outdir,'GHRSST_',predDate,'.grd',sep=''),overwrite=TRUE)
	writeRaster(tR2,filename=paste(outdir,'TurtleHabitat_',predDate,'.grd',sep=''),overwrite=TRUE)
}

## test the function for one date (predDate)
predDate <- dates[1]
outdir <- './'
varname <- 'SST'
min.var <- 17.5
max.var <- 18.5

test <- plot.turtle(predDate=predDate,outdir=outdir,varname=varname,min.var=min.var,max.var=max.var)

## apply the function over the list of dates using lapply 
## this will generate the plots into the specified directory

lapply(dates,FUN=plot.turtle,outdir=outdir,varname=varname,min.var=min.var,max.var=max.var)











