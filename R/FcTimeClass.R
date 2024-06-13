setClass("FcTime",representation(base="POSIXct",step="integer"),contains="Time",
	validity=function(object)
{
	if (any(is.na(object@base))) return("base NA")
	if (any(is.infinite(object@base))) return("base Inf")
	if (any(object@base < 0)) return("base < 0")

	if (any(is.na(object@step))) return("step NA")
	if (any(is.infinite(object@step))) return("step Inf")

	if (length(object@base) > 1 && length(object@step) > 1) {
		return("multiple base and step")
	}

	return(TRUE)
}
)
