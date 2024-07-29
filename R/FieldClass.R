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

sectioncs = function(field,long=c(0,360),lat=c(-90,90))
{
	sectioncsgrid(field,field@grid,long,lat)
}

zonalmean = function(field,mc.cores=1)
{
	zonalmeangrid(field@grid,field,mc.cores)
}

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
	} else if (method == "quad") {
		data = matrix(nrow=dim(field)[1],ncol=length(eta))
		dimnames(data)[[2]] = eta

		for (i in seq(along=eta)) {
			# choose levels for quadratic interp (lrr or llr)
			if (ind[i] == length(field@eta)-1) {
				ie = ind[i]+(-1):1
			} else {
				ie = ind[i]+0:2
				if (ind[i] > 1 && diff(range(field@eta[ie])) > diff(range(field@eta[ie-1]))) {
					ie = ie-1
				}
			}

			# Newton's polynom on non regular nodes
			data[,i] = newtonInterpv(field@eta[ie],field[,ie],eta[i])
		}
	} else {
		stop("method unsupported")
	}

	field@eta = eta
	setDataPart(field,data)
}
