#' @export
#' @import methods
#' @importClassesFrom gp.sa.diff DiffStatFrame
setClass("DBAStatFrame", contains="DiffStatFrame")

#' @export
#' @importClassesFrom gp.sa.diff DiffStatFrame
setClass("DGAStatFrame", contains="DiffStatFrame")
