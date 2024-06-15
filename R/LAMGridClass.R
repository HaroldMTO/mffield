setClass("LAMGrid",representation(nlong="integer",center="numeric"),contains="Grid",
	validity=function(object)
{
	nlat = length(object@nlong)
	if (nlat < 2) return("length(nlong) < 2")

	if (any(is.na(object@nlong))) return("nlong NA")
	if (any(is.infinite(object@nlong))) return("nlong Inf")
	if (any(object@nlong <= 0)) return("nlong <= 0")

	if (length(object@center) != 2) return("length(center) != 2")
	if (any(is.na(object@center))) return("center NA")
	if (any(is.infinite(object@center))) return("center Inf")
	if (abs(object@center[1]) > 90) return("center[1] (ie lat) out of [-90,90]")
	if ((object@center[2]+180)%%360 != object@center[2]+180) {
		return("center[2] (ie long) out of [-180,180[")
	}

	return(TRUE)
}
)

setMethod("==",signature(e1="LAMGrid",e2="LAMGrid"),def=function(e1,e2)
{
	if (! identical(e1@nlong,e2@nlong)) return(FALSE)

	if (! isTRUE(all.equal(e1@center,e2@center))) return(FALSE)

	callNextMethod(e1,e2)
}
)

setMethod("degrade","LAMGrid",def=function(grid,ilat,nlat,nlon)
{
	stopifnot(all(grid@nlong%%nlon == 0))
	nlatg = length(grid@nlong)

	if (missing(ilat)) {
		stopifnot(nlatg%%2 == 0)

		ilat = seq(1,nlatg,nlat)
	}

	stopifnot(all(ilat %in% seq(nlatg)) && identical(ilat,sort(ilat)))

	if (all(duplicated(grid@nlong)[-1])) {
		ndlon = grid@nlong[1]

		nlatg = length(grid@nlong)
		ilon = seq(1,ndlon,nlon)
		ind = rep((ilat-1)*ndlon,each=length(ilon))+ilon
		grid@lat = grid@lat[ind]
		grid@long = grid@long[ind]
		grid@nlong = rep(length(ilon),length(ilat))
	} else {
		# save
		nlong = grid@nlong

		clats = c(0,cumsum(grid@nlong))
		lats = grid@lat
		longs = grid@long

		# update
		nlon = as.integer(nlon)
		grid@nlong = grid@nlong[ilat]%/%nlon
		npdg = sum(grid@nlong)

		grid@lat = grid@long = numeric(npdg)

		off = 0
		for (i in ilat) {
			ind = clats[ilat[i]]+seq(1,nlong[i],nlon)
			np = length(ind)
			stopifnot(np == grid@nlong[i])

			grid@lat[off+1:np] = lats[ind]
			grid@long[off+1:np] = longs[ind]
			off = off+np
		}
	}

	grid
}
)

setMethod("extendPoints","LAMGrid",def=function(grid,x)
{
	x
}
)
