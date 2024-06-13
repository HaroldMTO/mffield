setClass("Grid",representation(lat="numeric",long="numeric",ind="integer"),
	validity=function(object)
{
	if (length(object@lat) != length(object@long)) return("length(lat) != length(long)")

	if (any(is.na(object@lat))) return("lat NA")
	if (any(is.infinite(object@lat))) return("lat Inf")
	if (any(abs(object@lat[1]) > 90)) return("lat out of [-90,90]")

	if (any(is.na(object@long))) return("long NA")
	if (any(is.infinite(object@long))) return("long Inf")
	if (any((object@long+180)%%360 != object@long+180)) return("long out of [-180,180[")

	#if (length(object@ind) > 0 && length(object@npdg0) == 0) return("ind with no npdg0")
	#if (length(object@ind) == 0 && length(object@npdg0) > 0) return("npdg0 with no ind")

	if (length(object@ind) > 0) {
		#if (length(object@npdg0) != 1) return("length(npdg0) != 1")
		#if (is.na(object@npdg0)) return("npdg0 NA")
		#if (is.infinite(object@npdg0)) return("npdg0 Inf")
		#if (object@npdg0 < 4) return("npdg0 < 4")

		if (any(is.na(object@ind))) return("ind NA")
		if (any(is.infinite(object@ind))) return("ind Inf")
		if (any(object@ind <= 0)) return("ind <= 0")
		#if (any(object@ind > object@npdg0)) return("ind > npdg0")
	}

	return(TRUE)
}
)

setGeneric("select",def=function(grid,npmax)
{
	standardGeneric("select")
}
)

setGeneric("degrade",signature="grid",def=function(grid,ilat,nlat=2,nlon=2)
{
	standardGeneric("degrade")
}
)

setGeneric("interpgrid",signature=c("g","grid"),def=function(g,grid,field,method,mc.cores)
{
	standardGeneric("interpgrid")
}
)

setGeneric("extendPoints",def=function(grid,x)
{
	standardGeneric("extendPoints")
}
)

setMethod("[","Grid",def=function(x,i,...)
{
	x@lat = x@lat[i]
	x@long = x@long[i]

	x
}
)

setMethod("==",signature(e1="Grid",e2="Grid"),def=function(e1,e2)
{
	isTRUE(all.equal(e1@lat,e2@lat)) & isTRUE(all.equal(e1@long,e2@long))
}
)

setMethod("length","Grid",def=function(x)
{
	length(x@lat)
}
)

setMethod("select","Grid",def=function(grid,npmax)
{
	round(seq(1,length(grid),length.out=npmax))
}
)

inDomain = function(grid,domain)
{
	ilat = domain@ylim[1] <= grid@lat & grid@lat <= domain@ylim[2]

	# case xlim=[0,0] (whole globe) is in alternative
	xlim = (domain@xlim+180)%%360-180
	if (diff(xlim) > 0) {
		ilon = xlim[1] <= grid@long & grid@long <= xlim[2]
	} else {
		ilon = grid@long <= xlim[2] | grid@long >= xlim[1]
	}

	ilat & ilon
}

toMatrix = function(grid,x)
{
	nlat = length(grid@nlong)
	clats = c(0,cumsum(grid@nlong))
	data = matrix(nrow=max(grid@nlong),ncol=nlat)

	for (ilat in seq(nlat)) {
		ip = clats[ilat]
		np = grid@nlong[ilat]
		data[1:np,ilat] = x[ip+1:np]
	}

	data
}

interp1d0 = function(datao,clato,ilato,ilono1,e1)
{
	ip1 = clato[ilato]+ilono1
	ip2 = clato[ilato]+ilono1%%(clato[ilato+1]-clato[ilato])+1
	stopifnot(all(clato[ilato] < ip2 & ip2 <= clato[ilato+1]))
	ip = ifelse(e1 <= .5,ip1,ip2)

	datao[ip,]
}

interp0 = function(datao,clato,ilato,ilono1,ilono2,e0,e1,e2)
{
	if (e0 <= .5) {
		interp1d0(datao,clato,ilato,ilono1,e1)
	} else {
		interp1d0(datao,clato,ilato+1,ilono2,e2)
	}
}

interp1d1 = function(datao,clato,ilato,ilono1,e1)
{
	ip1 = clato[ilato]+ilono1
	ip2 = clato[ilato]+ilono1%%(clato[ilato+1]-clato[ilato])+1
	stopifnot(all(clato[ilato] < ip2 & ip2 <= clato[ilato+1]))
	datao[ip1,]+e1*(datao[ip2,]-datao[ip1,])
}

interp1 = function(datao,clato,ilato,ilono1,ilono2,e0,e1,e2)
{
	d1 = interp1d1(datao,clato,ilato,ilono1,e1)
	d2 = interp1d1(datao,clato,ilato+1,ilono2,e2)

	d1+e0*(d2-d1)
}
