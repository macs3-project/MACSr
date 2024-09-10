#' @import basilisk
env_macs <- BasiliskEnvironment("env_macs", pkgname="MACSr",
                                packages = c("python=3.10"),
                                pip = c("macs3==3.0.2"))
