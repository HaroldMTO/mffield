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

setGeneric("zonalmeangrid",signature="grid",def=function(grid,field,mc.cores)
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
		my[ilat] = abs(yi-mean(yfn[indy]))

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
		if (length(ilat) == 0) {
			ilat = which.min(my)
		} else if (length(ilat) > 1) {
			cat("--> few lats crossing:",length(ilat),length(na.omit(x1)),"\n")
			ilat = ilat[which.min(my[ilat])]
			cat("--> lat crossing:",ilat,"\n")
		}

		ip = clats[ilat]+seq(grid@nlong[ilat])
		if (diff(x) < 0) {
			ind =  xf[ip] <= x[2] | xf[ip] >= x[1]
		} else {
			ind = x[1] <= xf[ip] & xf[ip] <= x[2]
		}

		x1 = xf[ip[ind]]
		data1 = field[ip[ind],]
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
		ii = which(x1[ind] <= x[2] | x1[ind] >= x[1])
		cat("--> long:",x1[ind],"\n")
	} else {
		ind = order(x1,na.last=NA)
		ii = which(x[1] <= x1[ind] & x1[ind] <= x[2])
	}

	ind = ind[ii]
	stopifnot(all(! is.na(x1[ind])))

	# not true if x1/data1 have been reset because of few lats
	if (length(dim(data1)) == 3) dim(data1) = c(2*nlat,dim(field)[2])
	stopifnot(all(! is.na(data1[ind,1])))

	if (length(long) == 1) {
		list(lats=x1[ind],data=data1[ind,])
	} else {
		list(longs=x1[ind],data=data1[ind,])
	}
}
)

mapdom = function(dom,grid,data,main=NULL,breaks="Sturges",palette="YlOrRd",pch=20,
	mar=c(2,2,3,5),mgp=c(2,1,0),cex=.6,ppi=72,quiet=FALSE,...)
{
	if (ppi > 144) stop("ppi > 144")

	# mar must be set before calling map AND passed to map
	# (because map resets mar internally and resets it on exit: this is then fake mar!)
	par(mar=mar,mgp=mgp)
	l = mapxy(dom,mar=mar,new=TRUE)
	box()

	h = hist(data,breaks,plot=FALSE)

	f = par("fin")
	npmax = as.integer(min(prod(f*ppi/(4*cex)),.Machine$integer.max))
	if (length(data) < npmax/100) {
		cat("--> very few points, magnify plotting symbol (x2)\n")
		cex = 2*cex
	} else if (length(data) > 1.2*npmax) {
		cex = max(.2,round(cex*sqrt(npmax/length(data)),3))
		npmax = as.integer(min(prod(f*ppi/(4*cex)),.Machine$integer.max))
	}

	if (length(data) > 1.2*npmax) {
		if (! quiet) {
			cat("--> reducing xy points from",length(data),"to",npmax,"and cex to",cex,"\n")
		}

		ind = select(grid,npmax)

		b2 = cut(data[ind],h$breaks)

		ilost = which(h$counts > 0 & table(b2) == 0)
		if (length(ilost) > 0) {
			b = cut(data,h$breaks)
			ind1 = which(b %in% levels(b)[ilost])
			ind = c(ind,ind1)
			stopifnot(all(! duplicated(ind)))
			if (! quiet) {
				cat("--> selecting back",length(ind1),"lost points in",length(ilost),
					"data bins\n")
			}
		}

		grid = grid[ind]
		data = data[ind]
	}

	if (cex < .2) cex = .2

	br = h$breaks

	ind = findInterval(data,br,rightmost.closed=TRUE)
	rev = regexpr("\\+$",palette) < 0
	cols = hcl.colors(length(br),sub("\\+$","",palette),rev=rev)

	tind = table(ind)

	p = .Last.projection()
	if (p$projection == "") {
		for (i in as.integer(names(sort(tind,decreasing=TRUE)))) {
			ii = which(ind == i)
			points(grid@long[ii],grid@lat[ii],col=cols[ind[ii]],pch=pch,cex=cex,...)
		}
	} else {
		for (i in as.integer(names(sort(tind,decreasing=TRUE)))) {
			ii = which(ind == i)
			mp = mapproject(grid@long[ii],grid@lat[ii])
			points(mp$x,mp$y,col=cols[ind],pch=pch,cex=cex,...)
		}
	}

	levels = sprintf("% .3g",br)
	if (any(duplicated(levels))) levels = sprintf("% .4g",br)
	maplegend(levels,col=cols)

	lines(l)
	title(main)
}

maplegend = function(breaks,col,...)
{
	u = par("usr")
	p = par("plt")

	width = (1-p[2])/6*diff(u[1:2])/diff(p[1:2])
	x = u[2]+width/3

	nl = length(breaks)
	height = diff(u[3:4])
	dy = height/(nl-1)
	y = u[3]
	ybas = y + dy*(seq(nl-1)-1)
	yhaut = ybas + dy

	rect(x,ybas,x+width,yhaut,col=col,border=NA,xpd=TRUE)

	op = par(las=2,yaxt="s")
	axis(4,c(ybas[1],yhaut),breaks,tick=FALSE,pos=x+width/12,mgp=c(1,.7,0),...)
	par(op)
}
