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

	if (length(object@ind) > 0) {
		if (any(is.na(object@ind))) return("ind NA")
		if (any(is.infinite(object@ind))) return("ind Inf")
		if (any(object@ind <= 0)) return("ind <= 0")
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

setGeneric("zonalmean",signature="g",def=function(g,field)
{
	standardGeneric("zonalmean")
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

setGeneric("sectiongeogrid",signature="grid",def=function(grid,field,long,lat)
{
	standardGeneric("sectiongeogrid")
}
)

setGeneric("sectioncsgrid",signature="grid",def=function(field,grid,long,lat)
{
	standardGeneric("sectioncsgrid")
}
)

setGeneric("zonalmeangrid",signature="grid",def=function(grid,field)
{
	standardGeneric("zonalmeangrid")
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

decircleLong = function(x)
{
	il = which(abs(diff(x)) > 180)
	stopifnot(length(il) <= 2)
	for (i in rev(il)) x[-(1:i)] = x[-(1:i)]-360*sign(diff(x)[i])

	x
}

decircleLat = function(x)
{
	il = which(abs(diff(x)) > 90)
	stopifnot(length(il) <= 2)
	for (i in rev(il)) x[-(1:i)] = x[-(1:i)]-180*sign(diff(x)[i])

	x
}

setMethod("sectiongeogrid","Grid",def=function(grid,field,long,lat)
{
	stopifnot(all(-90 <= lat & lat <= 90))

	# grid longs belong to [-180,180[
	if (length(long) == 1) {
		stopifnot(diff(lat) > 0)
		xf = grid@lat
		yf = grid@long
		x = lat
		y = (long+180)%%360-180
	} else if (length(lat) == 1) {
		stopifnot(diff(long) > 0)
		xf = grid@long
		yf = grid@lat
		x = (long+180)%%360-180
		y = lat
	} else {
		stop("lat or long must be of length 1")
	}

	nlat = length(grid@nlong)
	clats = c(0,cumsum(grid@nlong))

	data1 = array(dim=c(2,nlat,dim(field)[2]))
	x1 = array(dim=c(2,nlat))
	my = numeric(nlat)

	for (ilat in seq(nlat)) {
		off = clats[ilat]
		np = grid@nlong[ilat]
		ip = off+1:np
		xfn = xf[ip]
		yfn = yf[ip]
		yi = y
		# 1 more point (1st one) for global grid
		xfn = extendPoints(grid,xfn)
		yfn = extendPoints(grid,yfn)
		if (length(long) == 1) {
			xfn = decircleLat(xfn)
			yfn = decircleLong(yfn)
			stopifnot(all(diff(yfn) > -180))
			if (yi < min(yfn)) yi = yi+360
		} else {
			xfn = decircleLong(xfn)
			yfn = decircleLat(yfn)
			stopifnot(all(diff(xfn) > -180))
		}

		np = length(yfn)
		indy = which.min(abs(yi-yfn))+seq(-2,2)
		indy = indy[0 < indy & indy <= np]
		my[ilat] = mean(yfn[indy])

		# increasing or decreasing point values
		ind = yfn[-np] <= yi & yi < yfn[-1] | yfn[-1] <= yi & yi < yfn[-np]
		if (all(! ind)) next

		i1 = which(ind)
		stopifnot(length(i1) <= 2)
		stopifnot(all(i1 < np))

		# point np is only valid for global Gauss grids but never happens (i1 < np)
		i2 = i1+1
		e = (yi-yfn[i1])/(yfn[i2]-yfn[i1])
		# e can be 1 because of precision:
		# e = (b*(1-eps)-a)/(b-a)=1-b*eps/(b-a)=1-eps/(1-a/b), and e=1 for some cases
		stopifnot(all(0 <= e & e <= 1))

		for (i in seq(along=i1)) {
			data1[i,ilat,] = (1-e[i])*field[off+i1[i],]+e[i]*field[off+i2[i],]
			x1[i,ilat] = (1-e[i])*xfn[i1[i]]+e[i]*xfn[i2[i]]
		}
	}

	x1 = as.vector(x1)
	if (length(na.omit(x1)) < 5 && length(lat) == 1) {
		#stopifnot(length(lat) == 1)

		ilat = 1+(which(! is.na(x1))-1)%/%2
		cat("--> few lats crossing:",length(ilat),length(na.omit(x1)),"\n")
		if (length(ilat) == 0) {
			#ilat = which.min(abs(yi-my))
			ilat = which.min(my)
		} else if (length(ilat) > 1) {
			#ilat = ilat[which.min(abs(yi-my[ilat]))]
			ilat = ilat[which.min(my[ilat])]
		}

		cat("--> lat crossing:",ilat,"\n")
		ind = clats[ilat]+seq(grid@nlong[ilat])
		return(list(longs=xf[ind],lat=lat,data=field[ind,]))
	}

	if (all(is.na(x1))) {
		if (length(long) == 1) stop(paste("no lat/long crossing given long",long))
		if (length(lat) == 1) stop(paste("no lat/long crossing given lat",lat))
	}

	if (length(lat) == 1 && diff(x) < 0) {
		xn = x1[x1 < 0]
		xp = x1[x1 >= 0]
		indn = order(xn,na.last=NA)
		indp = order(xp,na.last=NA)
		ind = match(c(xp[indp],xn[indn]),x1)
		ii = which(x[1] <= x1[ind] | x1[ind] <= x[2])
		cat("--> long:",x1[ind],"\n")
	} else {
		ind = order(x1,na.last=NA)
		ii = which(x[1] <= x1[ind] & x1[ind] <= x[2])
	}

	ind = ind[ii]
	stopifnot(all(! is.na(x1[ind])))

	dim(data1) = c(2*nlat,dim(field)[2])
	stopifnot(all(! is.na(data1[ind,1])))

	if (length(long) == 1) {
		list(lats=x1[ind],data=data1[ind,])
	} else {
		list(longs=x1[ind],data=data1[ind,])
	}
}
)

