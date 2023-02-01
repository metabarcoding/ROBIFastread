#' @importFrom glue glue
#' @importFrom rlang is_true
NULL


#' Sets the verbose status
#'
#' @param state a logical value indicating the new verbose mode
#'
#' @export
#'
#' @examples
#'
#' state <- getOption("ROBITools2.verbose")
#' is_robi_verbose()
#' set_robi_verbose(TRUE)
#' is_robi_verbose()
#' set_robi_verbose(FALSE)
#' is_robi_verbose()
#' set_robi_verbose(state)
#'
set_robi_verbose <- function(state) {
  robiassert_arg(is.logical(state),
    arg = "state",
    message = "value must be a logical"
  )
  robiassert_arg(vec_size(state) == 1,
    arg = "state",
    message = "value must have a length of one"
  )

  options("ROBITools2.verbose" = state)
}

#' Indicates if ROBITools send messages to the terminal
#'
#' Mainly indicates if messages are sent to the user terminal
#' to indicate computation progression.
#'
#' @return a logical value
#' @export
#'
#' @rdname set_robi_verbose
#' @examples
#' is_robi_verbose()
#'
is_robi_verbose <- function(.envir = parent.frame()) {
  verbose <- tryCatch(get("verbose", envir = .envir),
    error = function(e) {
      FALSE
    }
  )

  (verbose ||
    rlang::is_true(getOption("ROBITools2.verbose"))) &&
    interactive() &&
    !isTRUE(getOption("rstudio.notebook.executing")) &&
    !isTRUE(getOption("knitr.in.progress"))
}


#' Writes a message on the console in verbose mode
#'
#' The message is written to the console if the code
#' is executed in the verbose mode. The verbose mode
#' can be specified by declaring a `verbose` logical
#' value to `TRUE`. If no `verbose` variable is defined,
#' the result of the `ROBITools2::is_robi_verbose`
#' function is used to decide of the verbose status.
#'
#'
#' @param message a string template sent to the `glue::glue` preprocessor.
#' @param ... other values to be concatenated to the end of the message
#' @param .envir the environment passed to the  `glue::glue` preprocessor
#'
#' @md
#' @export
#'
#' @examples
#' verbose <- TRUE
#' robimessage("Hello world !")
#' verbose <- FALSE
#' robimessage("Hello world !")
#' rm(verbose)
#' v_mode <- is_robi_verbose()
#' v_mode
#' robimessage("Hello world !")
#' set_robi_verbose(FALSE)
#' is_robi_verbose()
#' robimessage("Hello world !")
#' set_robi_verbose(TRUE)
#' is_robi_verbose()
#' robimessage("Hello world !")
#' set_robi_verbose(v_mode)
robimessage <- function(message, ..., .envir = parent.frame()) {
  if (is_robi_verbose(.envir)) {
    message(glue(message, ..., .envir = .envir))
  }
}
