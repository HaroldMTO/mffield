setClass("Domain",representation(xlim="numeric",ylim="numeric",proj="character",
	param="character"),
	validity=function(object)
{
	if (length(object@xlim) != 2) return("length(xlim) != 2")
	if (length(object@ylim) != 2) return("length(ylim) != 2")
	if (length(object@proj) > 1) return("length(proj) > 1")

	if (any(is.na(object@xlim))) return("xlim NA")
	if (any(is.na(object@ylim))) return("ylim NA")

	if (any(is.infinite(object@xlim))) return("xlim Inf")
	if (any(is.infinite(object@ylim))) return("ylim Inf")

	i = which(object@xlim != 180)
	if (any((object@xlim[i]+180)%%360 != object@xlim[i]+180))
		return("xlim out of [-180,180]")
	if (any(abs(object@ylim) > 90)) return("ylim out of [-90,90]")
	if (diff(object@ylim) <= 0) return("ylim not strictly increasing")

	return(TRUE)
}
)

setMethod("initialize","Domain",def=function(.Object,long="numeric",lat="numeric",...)
{
	if (! missing(long)) .Object@xlim = rangeLong(long)
	if (! missing(lat)) .Object@ylim = range(lat)

	.Object = callNextMethod(.Object,...)
	validObject(.Object)
	.Object
}
)

rangeLong = function(long)
{
	if (length(long) == 0) stop("no long provided")
	if (length(long) == 1) return(long)

	if (length(long) == 2) {
		dlong1 = long%%360
		dlong2 = (long+180)%%360-180
		if (diff(dlong1) < 0) {
			return(dlong2)
		} else {
			return(dlong1)
		}
	}

	#dlong1 = range(long%%360)-180
	dlong1 = range(long%%360)
	dlong2 = range((long+180)%%360-180)

	if (FALSE && diff(dlong1) > 180 && diff(dlong2) > 180) {
		if (diff(dlong1) > diff(dlong2)) {
			return(dlong1)
		} else {
			return(dlong2)
		}
	}

	if (diff(dlong1) <= diff(dlong2)) {
		return(dlong1)
	} else {
		return(dlong2)
	}
}

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

area = function(dom1,dom2)
{
	diff(dom1@xlim)*diff(dom1@ylim)/(diff(dom2@xlim)*diff(dom2@ylim))
}
