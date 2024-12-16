# divided difference of order n (recursive computation)
ddiffn = function(x,y,n,yprev)
{
	stopifnot(is.matrix(y))

	if (n == 0) {
		return(y)
	} else if (n == 1) {
		return((y[,2]-y[,1])/diff(x[1:2]))
	}

	if (missing(yprev)) {
		yx1 = ddiffn(x[-(n+1)],y[,-(n+1),drop=F],n-1)
	} else {
		yx1 = yprev
	}

	yx2 = ddiffn(x[-1],y[,-1,drop=F],n-1)
	(yx2-yx1)/(x[n+1]-x[1])
}

# Newton's polynom (or Newton's form of the Lagrange polynom)
Pn = function(x,y,xh)
{
	if (! is.matrix(y)) y = t(y)

	yprev = y
	for (i in seq(along=x)[-1]) yprev[,i] = ddiffn(x[1:i],y[,1:i,drop=F],i-1,yprev[,i-1])

	yprev = yprev[,-1]
	if (missing(xh)) return(yprev)

	y[,1]+sum(sapply(seq(along=x)[-length(x)],function(i) yprev[,i]*prod(xh-x[1:i])))
}

# Newton's polynom (or Newton's form of the Lagrange polynom)
Pnv = function(x,y,xh)
{
	if (! is.matrix(y)) y = t(y)
	nx = length(x)
	stopifnot(dim(y)[2] == nx)

	yprev = y
	for (i in seq(nx)[-1]) yprev[,i] = ddiffn(x[1:i],y[,1:i,drop=F],i-1,yprev[,i-1])

	yprev = yprev[,-1]
	if (missing(xh)) return(yprev)

	data = matrix(nrow=dim(y)[1],ncol=length(xh))
	for (j in seq(along=xh)) {
		data[,j] = y[,1]+sum(sapply(seq(nx-1),function(i) yprev[,i]*prod(xh[j]-x[1:i])))
	}

	data
}

newtonInterpv = function(x,y,xs=c())
{
	stopifnot(is.numeric(x),is.numeric(y))
	if (length(xs) != 0 && ! is.numeric(xs)) {
		stop("Argument 'xs' must be empty or a numeric vector.")
	}

	n = length(x)
	if (! is.matrix(y)) y = t(y)
	if (dim(y)[2] != n) stop("Vectors 'x' and 'y' must be of the same length.")

	# recursive computation of divided differences
	p = y
	for (k in 2:n) {
		for (l in k:n) p[,l] = (p[,l]-p[,k-1])/(x[l]-x[k-1])
	}

	if (length(xs) == 0) return(p)

	ys = matrix(nrow=dim(y)[1],ncol=length(xs))
	ys[,] = p[,n]

	# recursive (Horner's method) evaluation of Newton's polynom at xs
	for (k in 1:(n-1)) {
		for (j in seq(xs)) ys[,j] = p[,n-k]+(xs[j]-x[n-k])*ys[,j]
	}

	ys
}

dprodfreg = function(x,xh)
{
	if (length(x) == 1) return(1)

	# sum of Lagrange polynoms (li = prod(x-xj) with j!=i), really?
	sum(sapply(seq(along=x),function(i) prod(xh-x[-i])))
}

dprodf = function(x,y,xh)
{
	if (length(x) == 1) return(y)

	if (! is.matrix(y)) y = t(y)
	# sum of Lagrange polynoms (li = prod(x-xj) with j!=i), really?
	rowSums(sapply(seq(along=x),function(i) y[,i]*prod((xh-x[-i])/(x[i]-x[-i]))))
}

dPn = function(x,y,n,xh)
{
	newt = Pn(x,y,n)

	# derivative of Newton's polynom?
	sum(sapply(1:n,function(i) newt[i]*dprodfreg(x[1:i],xh)))
}
