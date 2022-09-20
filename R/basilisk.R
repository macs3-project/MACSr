#' @import basilisk
env_macs <- BasiliskEnvironment("env_macs", pkgname="MACSr",
                                packages = c("python=3.8"),
                                pip = c("macs3==3.0.0a7"))
