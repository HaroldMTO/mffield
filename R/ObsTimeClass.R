setClass("ObsTime",representation(date="POSIXct"),contains="Time",
	validity=function(object)
{
	if (any(is.na(object@date))) return("date NA")
	if (any(is.infinite(object@date))) return("date Inf")
	if (any(object@date < 0)) return("date < 0")

	return(TRUE)
}
)
