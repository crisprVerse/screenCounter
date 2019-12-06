#' @export
#' @import methods
#' @importClassesFrom gp.sa.diff DiffStatFrame
setClass("ScreenStatFrame", contains="DiffStatFrame")

#' @export
setClass("ScreenBarcodeStatFrame", contains="ScreenStatFrame")

#' @export
setClass("ScreenFeatureStatFrame", contains="ScreenStatFrame")
