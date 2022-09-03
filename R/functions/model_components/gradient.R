# logistic function to generate latitudinal temperature gradients based on four parameters
# accepts a vector of length 4 or a dataframe with four columns
gradient <- function(x, coeff, sdy = 0) { # parametrised with difference between cold and hot
  if(is.list(coeff) & !(is.data.frame(coeff) | is.matrix(coeff))) coeff = unlist(coeff)
  if(is.data.frame(coeff) | is.matrix(coeff)) {
    A = coeff[,1]
    DKA = coeff[,2]
    M = coeff[,3]
    B = coeff[,4]
    
    lat = t(data.frame(lat=x))
    lat = lat[rep(1, each=length(A)),]
    
    if(sdy == 0) {out = A + DKA/((1+(exp(B*(lat-M)))))
    } else {
      out = A + DKA/((1+(exp(B*(lat-M)))))+ rnorm(length(x),0,sdy)
    }
    
  } else {
    A = coeff[1]
    DKA = coeff[2]
    M = coeff[3]
    B = coeff[4]
    
    if(sdy == 0) {return(A + DKA/((1+(exp(B*(x-M))))))
    } else {
      out = A + DKA/((1+(exp(B*(x-M)))))+ rnorm(length(x),0,sdy)
    }
  }
  return(out)
}