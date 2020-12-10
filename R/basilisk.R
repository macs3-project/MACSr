#' @import basilisk
env_macs <- BasiliskEnvironment("env_macs", pkgname="MACSr",
                                packages = c("numpy>=1.17", "Cython>=0.29"),
                                pip = c("macs3"))
