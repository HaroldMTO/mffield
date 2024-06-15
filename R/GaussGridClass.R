setClass("GaussGrid",representation(nlong="integer",pole="numeric",stretch="numeric",
	theta="numeric"),contains="Grid",validity=function(object)
{
	nlat = length(object@nlong)
	if (nlat < 2) return("length(nlong) < 2")

	if (any(is.na(object@nlong))) return("nlong NA")
	if (any(is.infinite(object@nlong))) return("nlong Inf")
	if (any(object@nlong <= 0)) return("nlong <= 0")
	if (any(object@nlong[seq(nlat/2)] != rev(object@nlong[-seq(nlat/2)]))) {
		return("nlong not symmetric")
	}

	if (any(diff(object@nlong[seq(nlat/2)]) < 0)) return("nlong not increasing")

	if (length(object@pole) != 2) return("length(pole) != 2")
	if (any(is.na(object@pole))) return("pole NA")
	if (any(is.infinite(object@pole))) return("pole Inf")
	if (abs(object@pole[1]) > 90) return("pole[1] (ie lat) out of [-90,90]")
	if ((object@pole[2]+180)%%360 != object@pole[2]+180) {
		return("pole[2] (ie long) out of [-180,180[")
	}

	if (length(object@stretch) != 1) return("length(stretch) != 1")
	if (is.na(object@stretch)) return("stretch NA")
	if (is.infinite(object@stretch)) return("stretch Inf")
	if (object@stretch < 1) return("stretch < 1")

	if (length(object@theta) != nlat) return("length(theta) != nlat")
	if (any(is.na(object@theta))) return("theta NA")
	if (any(is.infinite(object@theta))) return("theta Inf")
	if (any(abs(object@theta) >= 90)) return("theta out of ]-90,90[")
	if (any(duplicated(object@theta)[-1])) return("theta duplicated")

	dlat = diff(object@theta)
	if (any(dlat > 0)) return("lat not strictly decreasing")
	#if (dlat[1] > 0 && any(dlat < 0)) return("lat not strictly increasing")

	return(TRUE)
}
)

setMethod("==",signature(e1="GaussGrid",e2="GaussGrid"),def=function(e1,e2)
{
	if (! identical(e1@nlong,e2@nlong)) return(FALSE)

	if (! isTRUE(all.equal(e1@pole,e2@pole))) return(FALSE)
	if (! isTRUE(all.equal(e1@stretch,e2@stretch))) return(FALSE)

	callNextMethod(e1,e2)
}
)

setMethod("!=",signature(e1="GaussGrid",e2="GaussGrid"),def=function(e1,e2)
{
	! e1 == e2
}
)

setAs("GaussGrid","LAMGrid",def=function(from)
{
	dom = new("Domain",long=from@long,lat=from@lat)
	xm = mean(dom@xlim)
	ym = mean(dom@ylim)
	new("LAMGrid",nlong=from@nlong,center=c(ym,xm),long=from@long,lat=from@lat)
}
)

setMethod("extendPoints","GaussGrid",def=function(grid,x)
{
	# 1 more point (1st one) for global grid
	c(x,x[1])
}
)

stretch = function(lat,stretch)
{
	if (stretch == 1) return(lat)

	stopifnot(stretch > 0)

	C2 = stretch^2
	mu = sin(lat*pi/180)
	sinla = ((C2+1)*mu+C2-1)/((C2-1)*mu+C2+1)
	theta = asin(sinla)*180/pi
	#Clat = seq(C,1/C,length.out=nlat)

	#schmidt = Clat*tan(acos(mu)/2)
	#theta2 = asin(cos(2*atan(schmidt)))*180/pi
	stopifnot(all(theta > lat))

	theta
}

csLong = function(grid)
{
	nlat = length(grid@nlong)
	longs = matrix(nrow=max(grid@nlong),ncol=nlat)

	for (i in seq(nlat)) {
		np = grid@nlong[i]
		longs[seq(np),i] = 180/pi*csLongi(np)
	}

	longs
}

csLongi = function(nlong)
{
	2*pi*seq(0,nlong-1)/nlong
}

pole = function(grid)
{
	c(180/pi*asin(grid@pole[1]),grid@pole[2])
}

compass = function(grid)
{
   nlat = length(grid@nlong)
   mucen = grid@pole[1]
   locen = pi/180*grid@pole[2]
   C = grid@stretch
   c2 = C^2

   sqm2 = sqrt(1-mucen^2)

   clats = c(0,cumsum(grid@nlong))

   l = m = numeric(sum(grid@nlong))

   for (ilat in seq(nlat)) {
      mu = sin(pi/180*grid@theta[ilat])
      sqmu2 = sqrt(1-mu^2)
      a = 1/(c2+1+(c2-1)*mu)
      b = c2-1+(c2+1)*mu

      nlon = grid@nlong[ilat]
      cslon = csLongi(nlon)
      csx = cos(cslon)
      csy = sin(cslon)
      ca1 = a*(2*C*mucen*sqmu2-b*sqm2*csx)
      gemu = a*(2*C*sqm2*sqmu2*csx+b*mucen)
      rcoslat = 1/sqrt(1-gemu^2)
      l[clats[ilat]+1:nlon] = -sqm2*csy*rcoslat
      m[clats[ilat]+1:nlon] = ca1*rcoslat
   }

   data.frame(l=l,m=m)
}

rotx = function(x,y,p,q)
{
   x*p-y*q
}

roty = function(x,y,p,q)
{
   x*q+y*p
}

geoCoords = function(grid)
{
   nlat = length(grid@nlong)
   mucen = grid@pole[1]
   locen = pi/180*grid@pole[2]
   C = grid@stretch
   c2 = C^2

   sqm2 = sqrt(1-mucen^2)

   clats = c(0,cumsum(grid@nlong))

   xp = cos(locen)
   yp = sin(locen)

   npdg = sum(grid@nlong)
   gelat = gelam = numeric(npdg)

   for (ilat in seq(nlat)) {
      mu = sin(pi/180*grid@theta[ilat])
		stopifnot(abs(mu) < 1)
      sqmu2 = sqrt(1-mu^2)
      a = 1/(c2+1+(c2-1)*mu)
      b = c2-1+(c2+1)*mu
      ccos = 2*C*sqmu2
      bcos = b*sqm2

      nlon = grid@nlong[ilat]
      cslon = csLongi(nlon)
      csx = cos(cslon)
      csy = sin(cslon)
      gemu = a*(ccos*sqm2*csx+b*mucen)
		stopifnot(all(abs(gemu) <= 1))

      gelat[clats[ilat]+1:nlon] = 180/pi*asin(gemu)
      loloc = rotx(mucen*csx,csy,xp,yp)
      lolos = roty(mucen*csx,csy,xp,yp)
		d = a/sqrt(1-gemu^2)
      geclo = d*(bcos*xp-ccos*loloc)
      geslo = d*(bcos*yp-ccos*lolos)
      sc = sign(geclo)
      c1 = min(sc,0)
      s2 = 2*min(sign(geslo),0)
      gelam[clats[ilat]+1:nlon] = 180/pi*(sc*asin(geclo)-pi*(c1+s2*(1+c1)))

		# reset long to 0 on poles
		ind = which(abs(gemu) == 1)
		gelam[clats[ilat]+ind] = 0
   }

	if (any(gelam < 0 | gelam >= 360)) {
		cat("--> shift longs:",range(gelam),"\n")
		gelam = gelam%%360
	}

   data.frame(lat=gelat,long=gelam)
}

setMethod("select","GaussGrid",def=function(grid,npmax)
{
	if (length(grid) <= npmax) return(seq(length(grid)))

	# note: 1 is trivial
	for (nbin in 2:20) {
      np = npdg%/%nbin
      ind = unlist(sapply(1:nbin,function(k) (k-1)*np+seq(1,np,by=nbin+1-k)))
      if (length(ind) <= npmax) break
   }

   ind
}
)

setMethod("degrade","GaussGrid",def=function(grid,ilat,nlat,nlon)
{
	stopifnot(all(grid@nlong%%nlon == 0))
	nlatg = length(grid@nlong)

	if (missing(ilat)) {
		stopifnot(nlatg%%2 == 0)

		ilatn = seq(1,nlatg/2,nlat)
		ilat = c(ilatn,rev(nlatg-ilatn+1))
	}

	stopifnot(all(ilat %in% seq(nlatg)) && identical(ilat,sort(ilat)))

	# save current values
	clats = c(0,cumsum(grid@nlong))
	lats = grid@lat
	longs = grid@long
	nlong = grid@nlong

	# update values
	grid@nlong = as.integer(1+(nlong[ilat]-1)%/%nlon)
	npdg = sum(grid@nlong)
	grid@theta = grid@theta[ilat]
	stopifnot(length(grid@nlong) == length(ilat))

	grid@lat = grid@long = numeric(npdg)

	off = 0
	for (i in seq(along=grid@nlong)) {
		ip = clats[ilat[i]]+seq(1,nlong[ilat[i]],nlon)
		np = length(ip)
		stopifnot(np == grid@nlong[i])

		grid@lat[off+1:np] = lats[ip]
		grid@long[off+1:np] = longs[ip]
		off = off+np
	}

	grid
}
)

ig = function(g,grid,field,method)
{
	nlato = length(g@nlong)

	# only for similar (Gauss) grids nested inside
	stopifnot(all(g@pole == grid@pole))
	stopifnot(nlato >= length(grid@nlong))
	stopifnot(max(g@nlong) >= max(grid@nlong))

	if (method == "linear") {
		interp1dfun = interp1d1
		interpfun = interp1
	} else if (method == "ppp") {
		interp1dfun = interp1d0
		interpfun = interp0
	} else {
		stop("unsupported method")
	}

	theta2 = stretch(grid@theta,grid@stretch)
	thetao2 = stretch(g@theta,g@stretch)
	stopifnot(thetao2[1] >= theta2[1])
	stopifnot(min(thetao2) <= min(theta2))

	data = matrix(nrow=sum(grid@nlong),ncol=dim(field)[2])
	dimnames(data) = list(NULL,dimnames(field)[[2]])

	clats = c(0,cumsum(grid@nlong))
	clato = c(0,cumsum(g@nlong))
	dlon = 360/g@nlong

	for (ilat in seq(along=grid@nlong)) {
		if (theta2[ilat] > thetao2[1]) {
			ilato = 1
		} else {
			ilato = max(which(thetao2 >= theta2[ilat]))
		}

		stopifnot(0 < ilato && ilato <= nlato)

		longs = 180/pi*csLongi(grid@nlong[ilat])
		ilon = seq(grid@nlong[ilat])
		e1 = longs/dlon[ilato]
		ilono1 = floor(e1)+1
		stopifnot(all(0 < ilono1 & ilono1 <= g@nlong[ilato]))
		e1 = e1-(ilono1-1)
		stopifnot(all(0 <= e1 & e1 < 1))

		ip = clats[ilat]+ilon
		if (theta2[ilat] > thetao2[1] || ilato == nlato) {
			data[ip,] = interp1dfun(field,clato,ilato,ilono1,e1)
			next
		}

		e0 = (thetao2[ilato]-theta2[ilat])/(thetao2[ilato]-thetao2[ilato+1])
		stopifnot(0 <= e0 && e0 < 1)
		if (e0 == 0) {
			data[ip,] = interp1dfun(field,clato,ilato,ilono1,e1)
			next
		}

		e2 = longs/dlon[ilato+1]
		ilono2 = floor(e2)+1
		stopifnot(all(0 < ilono2 & ilono2 <= g@nlong[ilato+1]))
		e2 = e2-(ilono2-1)
		stopifnot(all(0 <= e2 & e2 < 1))

		data[ip,] = interpfun(field,clato,ilato,ilono1,ilono2,e0,e1,e2)
	}

	data
}

interplat = function(ilat,g,grid,cs,clats,clato,dlono,theta2,thetao2,interp1dfun,
	interpfun,field)
{
	nlato = length(g@nlong)
	data = array(dim=c(grid@nlong[ilat],dim(field)[2]))

	for (i in seq(grid@nlong[ilat])) {
		ip = clats[ilat]+i
		if (cs$lat[ip] >= thetao2[1]) {
			ilato = 1
		} else {
			ilato = max(which(thetao2 >= cs$lat[ip]))
		}

		e1 = cs$long[ip]/dlono[ilato]
		ilono1 = floor(e1)+1
		stopifnot(0 < ilono1 && ilono1 <= g@nlong[ilato])
		e1 = e1-(ilono1-1)
		stopifnot(0 <= e1 && e1 < 1)

		# North from 1st lat or south from last lat (ie ilat = nlat)
		if (cs$lat[ip] >= thetao2[1] || ilato == nlato) {
			data[ip,] = interp1dfun(field,clato,ilato,ilono1,e1)
			next
		}

		e0 = (thetao2[ilato]-cs$lat[ip])/(thetao2[ilato]-thetao2[ilato+1])
		stopifnot(0 <= e0 && e0 < 1)
		if (e0 == 0) {
			data[ip,] = interp1dfun(field,clato,ilato,ilono1,e1)
			next
		}

		e2 = cs$long[ip]/dlono[ilato+1]
		ilono2 = floor(e2)+1
		stopifnot(0 < ilono2 && ilono2 <= g@nlong[ilato+1])
		e2 = e2-(ilono2-1)
		stopifnot(0 <= e2 && e2 < 1)

		data[ip,] = interpfun(field,clato,ilato,ilono1,ilono2,e0,e1,e2)
	}

	data
}

igpoint = function(g,grid,field,method,mc.cores)
{
	library(parallel)

	# consider the cs grid as a geo grid
	stopifnot(g@pole[1] == 1 && g@pole[2] == 0)

	nlato = length(g@nlong)

	if (method == "linear") {
		interp1dfun = interp1d1
		interpfun = interp1
	} else if (method == "ppp") {
		interp1dfun = interp1d0
		interpfun = interp0
	} else {
		stop("unsupported method")
	}

	theta2 = stretch(grid@theta,grid@stretch)
	thetao2 = stretch(g@theta,g@stretch)

	# find the geo coords, considered as cs coords in a geo grid (nsttyp=1)
	cs = geoCoords(grid)

	data = matrix(nrow=sum(grid@nlong),ncol=dim(field)[2])
	dimnames(data) = list(NULL,dimnames(field)[[2]])

	clats = c(0,cumsum(grid@nlong))
	clato = c(0,cumsum(g@nlong))
	dlono = 360/g@nlong

	ldata = mclapply(seq(along=grid@nlong),interplat,g,grid,cs,clats,clato,dlono,
		theta2,thetao2,interp1dfun,interpfun,field,mc.preschedule=FALSE,mc.cores=mc.cores)
	#for (ilat in seq(along=grid@nlong)) {
	#	ip = clats[ilat]+seq(grid@nlong[ilat])
	#	data[ip,] = interplat(ilat,g,grid,cs,clats,clato,dlono,theta2,thetao2,field)
	#}

	data = array(dim=c(sum(grid@nlong),dim(field)[2]))
	for (ilat in seq(along=grid@nlong)) {
		ip = clats[ilat]+seq(grid@nlong[ilat])
		data[ip,] = ldata[[ilat]]
	}

	data
}

igeopoint = function(g,grid,field,method)
{
	# consider the cs grid as a geo grid
	stopifnot(grid@pole[1] == 1 && grid@pole[2] == 0)

	nlat = length(grid@nlong)

	theta2 = stretch(grid@theta,grid@stretch)
	stopifnot(all(diff(theta2) < 0))
	thetao2 = stretch(g@theta,g@stretch)
	stopifnot(all(diff(thetao2) < 0))

	# find the geo coords, considered as cs coords in a geo grid (nsttyp=1)
	csg = geoCoords(g)

	data = matrix(0,nrow=sum(grid@nlong),ncol=dim(field)[2])
	dimnames(data) = list(NULL,dimnames(field)[[2]])
	wp = numeric(sum(grid@nlong))

	clats = c(0,cumsum(grid@nlong))
	clato = c(0,cumsum(g@nlong))
	dlon = 360/grid@nlong

	for (ilato in seq(along=g@nlong)) {
		for (i in seq(g@nlong[ilato])) {
			ip = clato[ilato]+i
			if (csg$lat[ip] >= theta2[1]) {
				ilat = 1
			} else {
				ilat = max(which(theta2 >= csg$lat[ip]))
			}

			# North from 1st lat or south from last lat (ie ilat = nlat)
			if (csg$lat[ip] >= theta2[1] || ilat == length(grid@nlong)) {
				e0 = 0
			} else {
				e0 = (theta2[ilat]-csg$lat[ip])/(theta2[ilat]-theta2[ilat+1])
				stopifnot(0 <= e0 && e0 < 1)
			}

			e1 = csg$long[ip]/dlon[ilat]
			ilon1 = floor(e1)+1
			stopifnot(0 < ilon1 && ilon1 <= grid@nlong[ilat])
			e1 = e1-(ilon1-1)
			stopifnot(0 <= e1 && e1 < 1)

			j1 = clats[ilat]+ilon1
			j2 = clats[ilat]+ilon1%%grid@nlong[ilat]+1
			data[j1,] = data[j1,]+(1-e0)*(1-e1)*field[ip,]
			data[j2,] = data[j2,]+(1-e0)*e1*field[ip,]
			wp[j1] = wp[j1]+(1-e0)*(1-e1)
			wp[j2] = wp[j2]+(1-e0)*e1

			if (e0 == 0) next

			e1 = csg$long[ip]/dlon[ilat+1]
			ilon1 = floor(e1)+1
			stopifnot(0 < ilon1 && ilon1 <= grid@nlong[ilat+1])
			e1 = e1-(ilon1-1)
			stopifnot(0 <= e1 && e1 < 1)

			j1 = clats[ilat+1]+ilon1
			j2 = clats[ilat+1]+ilon1%%grid@nlong[ilat+1]+1
			data[j1,] = data[j1,]+e0*(1-e1)*field[ip,]
			data[j2,] = data[j2,]+e0*e1*field[ip,]
			wp[j1] = wp[j1]+e0*(1-e1)
			wp[j2] = wp[j2]+e0*e1
		}
	}

	ind = wp > 0
	data[ind] = data[ind]/wp[ind]
	data
}

setMethod("interpgrid",signature(g="GaussGrid",grid="GaussGrid"),
	def=function(g,grid,field,method,mc.cores)
{
	if (all(grid@pole == g@pole)) {
		data = ig(g,grid,field,method)
	} else if (all(g@pole == c(1,0))) {
		cat("--> interpolation from a geo. grid to a CS grid\n")
		data = igpoint(g,grid,field,method,mc.cores)
	} else if (all(grid@pole == c(1,0))) {
		cat("--> interpolation from a CS grid to a geo grid: geoproject\n")
		data = igeopoint(g,grid,field,method)
	} else {
		cat("--> interpolation from a CS grid to a CS grid: geoproject + CSinterpolate\n")
		# create the equivalent geo grid
		ggeo = grid
		ggeo@pole = c(1,0)
		ggeo@stretch = 1
		clats = c(0,cumsum(ggeo@nlong))

		for (ilat in seq(along=ggeo@nlong)) {
			ip = clats[ilat]+seq(ggeo@nlong[ilat])
			ggeo@lat[ip] = ggeo@theta[ilat]
			ggeo@long[ip] = 180/pi*csLongi(ggeo@nlong[ilat])
		}

		cat("--> project on the geo grid\n")
		datageo = igeopoint(g,ggeo,field,method)

		cat("--> interpolate to the CS grid\n")
		data = igpoint(ggeo,grid,datageo,method,mc.cores)
	}

	data
}
)

setMethod("zonalmeangrid",signature("GaussGrid"),def=function(grid,field)
{
	nlat = length(grid@nlong)
	data = matrix(nrow=nlat,ncol=dim(field)[2])
	dimnames(data)[[2]] = dimnames(f)[[2]]

	clats = c(0,cumsum(grid@nlong))

	if (grid@pole[1] == 1) {
		for (ilat in seq(nlat)) {
			ip = clats[ilat]+seq(grid@nlong[ilat])
			data[ilat,] = apply(field[ip,,drop=FALSE],2,mean)
		}
	} else {
		# among lats, find for the same geo theta as the unstretched cs theta
		geolat = c(sort(grid@theta),90)
		ind = findInterval(grid@lat,geolat,all.inside=TRUE)
		stopifnot(all(0 < ind & ind <= nlat))

		for (ilat in seq(nlat)) data[ilat,] = apply(field[ind==ilat,,drop=FALSE],2,mean)
	}

	data
}
)

# divided difference of order n (recursive computation)
ddiffn = function(x,y,n,yprev)
{
	if (n == 0) {
		return(y)
	} else if (n == 1) {
		return(diff(y[1:2])/diff(x[1:2]))
	}

	if (missing(yprev)) {
		yx1 = ddiffn(x[-(n+1)],y[-(n+1)],n-1)
	} else {
		yx1 = yprev
	}

	yx2 = ddiffn(x[-1],y[-1],n-1)
	(yx2-yx1)/(x[n+1]-x[1])
}

# Newton's polynom (or Newton's form of the Lagrange polynom)
Pn = function(x,y,n,xh)
{
	yprev = y
	for (i in (1:n)+1) yprev[i] = ddiffn(x[1:i],y[1:i],i-1,yprev[i-1])

	yprev = yprev[-1]
	if (missing(xh)) return(yprev)

	y[1]+sum(sapply(1:n,function(i) yprev[i]*prod(xh-x[1:i])))
}

dprodf = function(x,xh)
{
	if (length(x) == 1) return(1)

	# sum of Lagrange polynoms (li = prod(x-xj) with j!=), really?
	sum(sapply(seq(along=x),function(j) prod(xh-x[-j])))
}

dPn = function(x,y,n,xh)
{
	newt = Pn(x,y,n)

	# derivative of Newton's polynom?
	sum(sapply(1:n,function(i) newt[i]*dprodf(x[1:i],xh)))
}

gradl = function(data,grid,method="newton",mc.cores=1)
{
	stopifnot(is(grid,"GaussGrid"))
	stopifnot(method %in% c("newton","FD"))

	clats = c(0,cumsum(grid@nlong))
	nlat = length(grid@nlong)

	datal = data

	for (i in seq(nlat)) {
		ip = clats[i]+seq(grid@nlong[i])
		np = length(ip)

		# circular points (1 or 2 shifted)
		ip1 = c(ip[np],ip[-np])
		ip3 = c(ip[-1],ip[1])
		ip0 = c(ip[np-1:0],ip[0:1-np])
		ip4 = c(ip[-(1:2)],ip[1:2])

		if (method == "newton") {
			m = t(matrix(c(ip0,ip1,ip,ip3,ip4),ncol=5))

			for (j in seq(dim(data)[2])) {
				a = mclapply(seq(dim(m)[2]),function(k) dPn(1:5,data[m[,k],j],4,3),
					mc.cores=mc.cores)
				datal[ip,j] = simplify2array(a)
			}
		} else {
			dlon = abs(diff(grid@long[ip[1:2]]))
			a = cos(grid@lat[ip[1]]*pi/180)
			dl = sqrt(diff(grid@lat[ip[1:2]])^2+(a*pmin(dlon,360-dlon))^2)
			dd = data[ip3,]-data[ip1,]
			#d3 = dd/(2*dl)
			#d5 = dd/(4*dl)
			#data = 2*d3-d5
			datal[ip,] = 3/4*dd/dl
		}
	}

	datal
}

setMethod("sectioncsgrid",signature(grid="GaussGrid"),
	def=function(field,grid,long=c(0,360),lat=c(-90,90))
{
	stopifnot(all(-90 <= lat & lat <= 90))

	clats = c(0,cumsum(grid@nlong))
	dlon = 360/grid@nlong

	if (length(long) == 1) {
		if (length(lat) == 1) {
			dlat = 180/(length(grid@nlong)+1)
			ilat = which.max(grid@theta <= lat)
			e0 = (90-dlat/2-grid@theta[ilat])/dlat-ilat+1
			ind = c(ilat,ilat+1)
		} else {
			ind = which(min(lat) <= grid@theta & grid@theta <= max(lat))
		}

		data = matrix(nrow=length(ind),ncol=dim(field)[2])
		lats = longs = numeric(length(ind))

		# grid (CS) longs belong to [0,360[ (use positive values)
		long = long%%360

		for (i in seq(along=ind)) {
			ilat = ind[i]
			e = long/dlon[ilat]
			ilon = floor(e)+1
			e = e-(ilon-1)
			ip1 = clats[ilat]+ilon
			ip2 = clats[ilat]+ilon%%(clats[ilat+1]-clats[ilat])+1
			data[i,] = (1-e)*field[ip1,]+e*field[ip2,]
			lats[i] = (1-e)*grid@lat[ip1]+e*grid@lat[ip2]
			longs[i] = grid@long[ip1]
		}

		if (length(lat) == 1) data = data[1,]+e0*(data[2,]-data[1,])

		list(theta=grid@theta[ind],lats=lats,longs=longs,data=data)
	} else if (length(lat) == 1) {
		# cannot sort or do modulo because long is used as ordered values
		long = unique(long)
		stopifnot(length(long) == 2)
		stopifnot(all(0 <= long & long <= 360))

		ilat = max(which(grid@theta >= lat))
		longs = 180/pi*csLongi(grid@nlong[ilat])

		if (long[1] < long[2]) {
			ind = which(long[1] <= longs & longs <= long[2])
		} else {
			ind = c(which(long[1] <= longs & longs <= 360),
				which(0 <= longs & longs <= long[2]))
		}

		data = field[clats[ilat]+ind,j]

		#e = (lat-grid@theta[ilat])/diff(grid@theta[ilat+0:1])

		ip = clats[ilat]+ind

		list(lambda=longs[ind],lats=grid@lat[ip],longs=grid@long[ip],data=data)
	}
}
)

