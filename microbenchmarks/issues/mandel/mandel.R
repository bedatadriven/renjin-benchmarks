## mandel
## Discussed in http://stackoverflow.com/questions/9968578/speeding-up-julias-poorly-written-r-examples/9969530#9969530

mandel = function(z) {
    c = z
    maxiter = 80
    for (n in 1:maxiter) {
        if (Mod(z) > 2) return(n-1)
        z = z^2+c
    }
    return(maxiter)
}

mandelperf = function() {
    re = seq(-2,0.5,.1)
    im = seq(-1,1,.1)
    M = matrix(0.0,nrow=length(re),ncol=length(im))
    count = 1
    for (r in re) {
        for (i in im) {
            M[count] = mandel(complex(real=r,imag=i))
            count = count + 1
        }
    }
    return(M)
}

stopifnot(sum(mandelperf()) == 14791)

