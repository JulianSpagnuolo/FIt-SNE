fftRtsne <- function(X, 
		     dims=2, perplexity=30, theta=0.5,
		     check_duplicates=TRUE,
		     max_iter=1000,
		     fft_not_bh = TRUE,
		     ann_not_vptree = TRUE,
		     stop_lying_iter=250,
		     exaggeration_factor=12.0, no_momentum_during_exag=FALSE,
		     start_late_exag_iter=-1.0,late_exag_coeff=1.0,
		     n_trees=50, search_k = -1,rand_seed=-1,
		     nterms=3, intervals_per_integer=1, min_num_intervals=50, 
		     data_path=NULL, result_path=NULL,
		     fast_tsne_path=NULL, nthreads=NULL, ...)
{
  #' @export
  #' @title ffRtsne
  #' @author George Linderman, Cristoph H, Julian Spagnuolo
  #' @param X numeric matrix or dataframe containing data to be tsne-fied
  #' @param dims integer. number of output dimensions to reduce to. Default is 2.
  #' @param perplexity integer. Default is 30.
  #' @param theta numeric. Set to 0 for exact.  If non-zero, then will use either Barnes Hut or FIt-SNE based on fft_not_bh. If Barnes Hut is used, then this determines the accuracy of BH approximation. Default is 0.5
  #' @param max_iter integer. Maximum iterations to run tsne over the data.
  #' @param fft_not_bh logical. Whether to run the Fast Fourier Transform tSNE (TRUE), or if set to FALSE, will use the BH implementation (slower). Default is TRUE.
  #' @param ann_not_vptree logical. If TRUE will approximate nearest neighbors using [Annoy](https://github.com/spotify/annoy), if FALSE the original BH-tSNE method using vantage-point trees will be used. Default is TRUE.
  #' @param stop_lying_iter integer. Determines the iteration at which early exageration is switch off. Default is 250.
  #' @param exaggeration_factor numeric. Early exageration coefficient, altering can lead to improved embedding of swissrolls and other synthetic datasets. Default is 12.0
  #' @param no_momentum_during_exag logical. If FALSE, momentum and other optimizations will be used in the gradient decent phase. If TRUE, standard gradient descent will be used. Useful for testing large exaggeration coefficients. Default is FALSE.
  #' @param start_late_exag_iter numeric. Iteration at which to start late exageration, if set to -1, late exageration will not be used. Default is -1.0
  #' @param late_exag_coeff numeric. Increasing this parameter may improve separation of the clusters in the reduced space. Default is 1.0
  #' @param n_trees integer. Default is 50
  #' @param search_k numeric. Default is -1.
  #' @param rand_seed numeric. Set a seed number for reproducibility, if -1 the random seed will be set using the system time. Default is -1.
  #' @param intervals_per_integer integer. See Details, must be >0. Default is 1
  #' @param min_num_intervals integer. See Details. Default is 50
  #' @param data_path character. File path of a data file. If NULL, X will be written to a temp file. Default is NULL
  #' @param result_path character. File path to write results to, if NULL results will be written to a temp file prior to being read into R. Default is NULL
  #' @param fast_tsne_path character. File path to the location of the compiled fast_tsne_binary. This version of the forked git repo will compile the fast_tsne binary to the package installation directory, if this param is NULL then the function will default to this location. You can change this behaviour if you have a custom install of the original repo. Default is NULL
  #' @param nthreads integer. Number of threads used by computeGaussianPerplexity, if NULL, will use detectCores()-1 to determine the number of threads to use.
  #' @details 
  #' ann_not_vptree
  #'  Using "near" neighbors as opposed to strictly "nearest" neighbors is faster, but also has a smoothing effect,
  #'  which can be useful for embedding some datasets (see Linderman et al. (2017)).
  #'  If subtle detail is required (e.g. in identifying small clusters), then use vantage-point trees
  #'  (which is also multithreaded in this implementation).
  #' 
  #' min_num_intervals & intervals_per_integer
  #' 
  #' Let:
  #'    maxloc = ceiling(max(max(X)))
  #'    and
  #'    minloc = floor(min(min(X)))
  #'    i.e. the data points are in a minloc^n.dims by maxloc^n.dims interval/square. The number os intervals in each dimension
  #'    is either min_num_intervals or ceiling(maxloc - minloc)/intervals_per_integer whichever is larger.
  #'    
  #' @return Matrix of the reduced tSNE dimensions, ncol = dims, nrow = nrow(X).
  #' @references George C. Linderman, Manas Rachh, Jeremy G. Hoskins, Stefan Steinerberger, Yuval Kluger. (2017). Efficient Algorithms for t-distributed Stochastic Neighborhood Embedding. (2017) arXiv:1712.09005
  #' 
  #' @examples 
  #' # import the iris data set and remove any duplicates 
  #' data(iris)
  #' 
  #' u.iris <- unique(iris)
  #' 
  #' # For a standard BH-tsne
  #' 
  #' bh <- fftRtsne(X=u.iris[,1:4], dims = 2, fft_not_bh = F, max_iter = 1000, rand_seed=42)
  #' u.iris$bh.tsne1 <- bh[,1]
  #' u.iris$bh.tsne2 <- bh[,2]
  #' 
  #' ggplot(u.iris, aes(x=bh.tsne1, y=bh.tsne2, colour=Species)) +geom_point() +scale_color_colorblind() +theme_bw() +theme(aspect.ratio=1)
  #' 
  #' # For the fast Fourier transform tSNE
  #' 
  #' ff <- fftRtsne(X=u.iris[,1:4], dims=2, fft_not_bh=T, max_iter=1000, rand_seed=42)
  #' u.iris$ff.tsne1 <- ff[,1]
  #' u.iris$ff.tsne2 <- ff[,2]
  #' 
  #' ggplot(u.iris, aes(x=ff.tsne1, y=ff.tsne2, colour=Species)) +geom_point() +scale_color_colorblind() +theme_bw() +theme(aspect.ratio=1)
  #' 
  #' @importFrom parallel detectCores
	
  if (is.null(data_path)) {
		data_path <- tempfile(pattern='fftRtsne_data_', fileext='.dat')
	}
	if (is.null(result_path)) {
		result_path <- tempfile(pattern='fftRtsne_result_', fileext='.dat')
	}
	#if (is.null(fast_tsne_path)) { ## No longer necessary with default install location using makefile. JS
	#	fast_tsne_path <- system2('which', 'fast_tsne', stdout=TRUE)
	#}
  if(is.null(fast_tsne_path)) ### Added in to default to R package install location. JS.
  {
    fast_tsne_path <- paste(normalizePath(find.package("fftRtsne")), "fast_tsne", sep="/")
  }
	fast_tsne_path <- normalizePath(fast_tsne_path)
	if (!file_test('-x', fast_tsne_path)) {
		stop(fast_tsne_path, " does not exist or is not executable; check your fast_tsne_path parameter")
	}

	if(is.null(nthreads)) ## added in to make maximal use of detected cores, minus one for safety :) JS
	{
	  nthreads <- detectCores()-1
	}
	
	is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol

	if (!is.numeric(theta) || (theta<0.0) || (theta>1.0) ) { stop("Incorrect theta.")}
	if (nrow(X) - 1 < 3 * perplexity) { stop("Perplexity is too large.")}
	if (!is.matrix(X)) ## changed to convert to matrix if it is a dataframe, added verbosity if you tried to give it something else. JS
	  { 
	    if(is.data.frame(X))
	    {
	      X <- as.matrix(X)
	    }
	    else(stop("Not a data.frame or matrix"))
	  } 
	if (!(max_iter>0)) { stop("Incorrect number of iterations.")}
	if (!is.wholenumber(stop_lying_iter) || stop_lying_iter<0) { stop("stop_lying_iter should be a positive integer")}
	if (!is.numeric(exaggeration_factor)) { stop("exaggeration_factor should be numeric")}
	if (!is.wholenumber(dims) || dims<=0) { stop("Incorrect dimensionality.")}
	if (search_k == -1) { search_k = n_trees*perplexity*3 }

	if (fft_not_bh){
	  nbody_algo = 2;
	}else{
	  nbody_algo = 1;
	}
	
	if (ann_not_vptree){
	  knn_algo = 1;
	}else{
	  knn_algo = 2;
	}
	tX = c(t(X))

	f <- file(data_path, "wb")
	n = nrow(X);
	D = ncol(X);
	writeBin(as.integer(n), f,size= 4)
	writeBin( as.integer(D),f,size= 4)
	writeBin( as.numeric(0.5), f,size= 8) #theta
	writeBin( as.numeric(perplexity), f,size= 8) #theta
	writeBin( as.integer(dims), f,size=4) #theta
	writeBin( as.integer(max_iter),f,size=4)
	writeBin( as.integer(stop_lying_iter),f,size=4)
	writeBin( as.integer(-1),f,size=4) #K
	writeBin( as.numeric(-30.0), f,size=8) #sigma
	writeBin( as.integer(nbody_algo), f,size=4)  #not barnes hut
	writeBin( as.integer(knn_algo), f,size=4) 
	writeBin( as.numeric(exaggeration_factor), f,size=8) #compexag
	writeBin( as.integer(no_momentum_during_exag), f,size=4) 
	writeBin( as.integer(n_trees), f,size=4) 
	writeBin( as.integer(search_k), f,size=4) 
	writeBin( as.integer(start_late_exag_iter), f,size=4) 
	writeBin( as.numeric(late_exag_coeff), f,size=8) 
	
	writeBin( as.integer(nterms), f,size=4) 
	writeBin( as.numeric(intervals_per_integer), f,size=8) 
	writeBin( as.integer(min_num_intervals), f,size=4) 
	tX = c(t(X))
	writeBin( tX, f) 
	writeBin( as.integer(rand_seed), f,size=4) 
	close(f) 

	flag= system2(command=fast_tsne_path, args=c(data_path, result_path, nthreads));
	if (flag != 0) {
		stop('tsne call failed');
	}
	f <- file(result_path, "rb")
	initialError <- readBin(f, integer(), n=1, size=8);
	n <- readBin(f, integer(), n=1, size=4);
	d <- readBin(f, integer(), n=1,size=4);
	Y <- readBin(f, numeric(), n=n*d);
	Yout <- t(matrix(Y, nrow=d));
	close(f)
	file.remove(data_path)
	file.remove(result_path)
	Yout
}





