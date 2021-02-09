#' @import basilisk
env_macs <- BasiliskEnvironment("env_macs", pkgname="MACSr",
                                packages = c("numpy>=1.17"),
                                pip = c("macs3==3.0.0a6"))
