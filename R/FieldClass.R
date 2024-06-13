setClass("Field",representation(grid="Grid",eta="numeric"),contains="array",
	validity=function(object)
{
	if (length(object@.Data) > 1) {
		npdg = length(object@grid)
		nlev = length(object@eta)
		if (length(object@.Data) != npdg*nlev) return("dim(data) != [grid,eta]")

		if (! is.null(dim(object@.Data))) {
			if (! npdg %in% dim(object@.Data)) return("dim(data) inconsistent with pdg")
			if (! nlev %in% dim(object@.Data)) return("dim(data) inconsistent with levels")
			if (diff(match(c(npdg,nlev),dim(object@.Data))) < 0) {
				return("data may have swapped dimensions")
			}
		}
	}

	return(TRUE)
}
)

zoom = function(field,domain)
{
	# full grid only
	stopifnot(length(field@grid@ind) == 0)

	ind = inDomain(field@grid,domain)
	if (all(ind)) return(field)

	g = field@grid[ind]
	if (diff(domain@xlim) < 360) g = as(g,"LAMGrid")
	g@ind = which(ind)
	field@grid = g

	setDataPart(field,field[ind,,drop=FALSE])
}

setGeneric("sectiongeogrid",signature="grid",def=function(grid,field,long,lat)
{
	standardGeneric("sectiongeogrid")
}
)

setGeneric("sectioncs",signature="grid",def=function(field,grid,long,lat)
{
	standardGeneric("sectioncs")
}
)

interp = function(field,grid,method="linear",mc.cores=1)
{
	data = interpgrid(field@grid,grid,field,method,mc.cores)
	field@grid = grid
	setDataPart(field,data)
}


nextPoint = function(ilon,offgp,nlon)
{
	ip = offgp+(ilon%%nlon)+1
	stopifnot(all(offgp < ip))
	ip
}

sectiongeo = function(field,long=c(0,360),lat=c(-90,90))
{
	sectiongeogrid(field@grid,field,long,lat)
}

setMethod("sectiongeogrid","Grid",def=function(grid,field,long,lat)
{
	stopifnot(all(-90 <= lat & lat <= 90))

	# frame longs belong to [-180,180[
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

setMethod("sectioncs",signature(grid="GaussGrid"),def=function(field,grid,long=c(0,360),
	lat=c(-90,90))
{
	stopifnot(all(-90 <= lat & lat <= 90))

	clats = c(0,cumsum(frame$nlong))
	dlon = 360/frame$nlong

	if (length(long) == 1) {
		if (length(lat) == 1) {
			dlat = 180/frame$nlat
			ilat = which.max(frame$theta <= lat)
			e0 = (90-dlat/2-frame$theta[ilat])/dlat-ilat+1
			ind = c(ilat,ilat+1)
		} else {
			ind = which(min(lat) <= frame$theta & frame$theta <= max(lat))
		}

		data1 = matrix(nrow=length(ind),ncol=dim(data)[2])
		lats = longs = numeric(length(ind))

		# frame longs belong to [-180,180[
		long = (long+180)%%360-180

		for (i in seq(along=ind)) {
			ilat = ind[i]
			e = long/dlon[ilat]
			ilon = floor(e)+1
			e = e-(ilon-1)
			ip1 = clats[ilat]+ilon
			ip2 = clats[ilat]+ilon%%(clats[ilat+1]-clats[ilat])+1
			data1[i,] = (1-e)*data[ip1,]+e*data[ip2,]
			lats[i] = (1-e)*frame$lat[ip1]+e*frame$lat[ip2]
			longs[i] = frame$long[ip1]
		}

		if (length(lat) == 1) data1 = data1[1,]+e0*(data1[2,]-data1[1,])

		list(lat=frame$theta[ind],lats=lats,longs=longs,data=data1)
	} else if (length(lat) == 1) {
		long = unique(long)
		stopifnot(length(long) == 2)
		stopifnot(all(0 <= long & long <= 360))

		ilat = max(which(frame$theta >= lat))
		longs = 180/pi*csLongi(frame$nlon[ilat])

		if (long[1] < long[2]) {
			ind = which(long[1] <= longs[,ilat] & longs[,ilat] <= long[2])
		} else {
			ind = c(which(long[1] <= longs[,ilat] & longs[,ilat] <= 360),
				which(0 <= longs[,ilat] & longs[,ilat] <= long[2]))
		}

		data1 = matrix(nrow=length(ind),ncol=dim(data)[2])

		for (j in seq(dim(data)[2])) {
			datag = datatoGauss(data[,j],frame)
			data1[,j] = datag[ind,ilat]
		}

		lats = longs = numeric(length(ind))

		e = (lat-frame$theta[ilat])/diff(frame$theta[ilat+0:1])

		ip = clats[ilat]+ind

		list(long=longs[ind,ilat],lats=frame$lat[ip],longs=frame$long[ip],data=data1)
	}
}
)

interpAB = function(field,eta,method="linear")
{
	ind = findInterval(eta,field@eta)
	stopifnot(all(0 < ind & ind < length(field@eta)))
	e = (eta-field@eta[ind])/(field@eta[ind+1]-field@eta[ind])
	stopifnot(all(0 <= e & e < 1))

	if (method == "ppp") {
		ind = ifelse(e <= .5,ind,ind+1)
		data = field[,ind,drop=FALSE]
	} else if (method == "linear") {
		data = matrix(nrow=dim(field)[1],ncol=length(eta))
		dimnames(data)[[2]] = eta

		for (i in seq(along=eta)) {
			d1 = field[,ind[i]]
			if (e[i] < 1.e-4) {
				data[,i] = d1
				next
			}

			d2 = field[,ind[i]+1]
			if (1-e[i] < 1.e-4) {
				data[,i] = d2
				next
			}

			data[,i] = d1+e[i]*(d2-d1)
		}
	} else {
		stop("method unsupported")
	}

	field@eta = eta
	setDataPart(field,data)
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
