###
## Transfer function for calcite from Kim and O'Neil (1997),
#  re-expressed in Leng and Marshall (2004):
#
calcite_temp <- function(d18O_calcite,d18O_sw) {
  13.8 - 4.58 * (d18O_calcite - d18O_sw) + 0.08 * (d18O_calcite - d18O_sw)^2
}

###
## Transfer function for phosphate from Puceat et al. (2010),
#  which version?
#
phosphate_temp <- function(d18O_phosphate,d18O_sw) {

}
