#' install MACS python package
#'
#' @rdname install_MACS
#' @param envname Virtual environment to install the MACS python
#'     module
#' @return install MACS python module
#' @importFrom reticulate virtualenv_create virtualenv_install
#'     virtualenv_list use_virtualenv
#' @export

install_MACS <- function(envname = "MACS") {
    stopifnot(is.character(envname) && length(envname) == 1)
    is_windows <- identical(.Platform$OS.type, "windows")
    is_osx <- Sys.info()["sysname"] == "Darwin"
    is_linux <- identical(tolower(Sys.info()[["sysname"]]), "linux")
    if (!is_windows && !is_osx && !is_linux) {
        stop(
            "Unable to install 'MACS' on this platform. ",
            "Binary installation is available for Windows, macOS, and Linux"
        )
    }

    if (!envname %in% virtualenv_list()) {
        pkgs <- readLines("python_requirements.txt")
        virtualenv_create(envname)
        virtualenv_install(envname, pkgs)
    }
    use_virtualenv(virtualenv=envname)

    invisible(.MACS())
}
