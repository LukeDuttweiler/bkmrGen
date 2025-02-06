Matrix_to_slam <- function(Mat){
  p <- rep(seq_along(diff(Mat@p)), diff(Mat@p))
  Mat <- slam::simple_triplet_matrix(i = (Mat@i)+1, j = p, v = Mat@x,
                                   nrow = Mat@Dim[1], ncol = Mat@Dim[2])
  return(Mat)
}
