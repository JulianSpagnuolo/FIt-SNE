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
  #' @param X numeric matrix or dataframe containing data to be tsne-fied
  #' @param dims integer. number of output dimensions to reduce to. Default is 2.
  #' @param perplexity integer. Default is 30.
  #' @param theta numeric. Default is 0.5
  #' @param max_iter integer. Maximum iterations to run tsne over the data.
  #' @param fft_not_bh logical. Whether to run the Fast Fourier Transform tSNE (TRUE), or if set to FALSE, will use the BH implementation (slower). Default is TRUE.
  #' @param ann_not_vptree logical. Default is TRUE.
  #' @param stop_lying_iter integer. Default is 250.
  #' @param exaggeration_factor numeric. Default is 12.0
  #' @param no_momentum_during_exag logical. Default is FALSE.
  #' @param start_late_exag_iter numeric. Default is -1.0
  #' @param late_exag_coeff numeric. Default is 1.0
  #' @param n_trees integer. Default is 50
  #' @param intervals_per_integer integer. Default is 1
  #' @param min_num_intervals integer. Default is 50
  #' @param data_path character. Default is NULL
  #' @param result_path character. Default is NULL
  #' @param fast_tsne_path character. File path to the location of the compiled fast_tsne_binary. This version of the forked git repo will compile the fast_tsne binary to the package installation directory, if this param is NULL then the function will default to this location. You can change this behaviour if you have a custom install of the original repo. Default is NULL
  #' @param nthreads integer. Number of threads to use in the FF-tSNE, if NULL, will use detectCores()-1 to determine the number of threads to use.
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





