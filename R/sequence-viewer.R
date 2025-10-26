#' Create a sequence viewer output element
#' @param id The ID of the output element
#' @param width Width of the container
#' @param height Height of the container
#' @export
sequenceOutput <- function(id, width = "100%", height = "auto") {
  htmltools::tagList(
    # Container div that will be found by JavaScript
    htmltools::div(
      id = id,
      class = "sequence-viewer-container",
      style = sprintf("width: %s; height: %s; min-height: 100px;", width, height)
    )
  )
}

#' Render a sequence with quality scores
#' @param expr Expression that returns a list with bases, qscores, crlStart, and crlEnd
#' @export
renderSequence <- function(expr) {
  # Convert the expression to a function
  func <- shiny::exprToFunction(expr, env = parent.frame())
  
  function(session, name, ...) {
    # Try to get data
    tryCatch({
      value <- func()
      
      # Basic validation
      if (is.null(value)) {
        warning("No sequence data returned")
        return(NULL)
      }
      
      if (!is.list(value)) {
        warning("Invalid sequence data format")
        return(NULL)
      }
      
      if (!all(c("bases", "qscores", "crlStart", "crlEnd") %in% names(value))) {
        warning("Missing required sequence data fields")
        return(NULL)
      }
      
      # Type conversion
      value$bases <- unlist(value$bases)  # Handle potential list input
      value$qscores <- as.numeric(value$qscores)
      value$crlStart <- as.numeric(value$crlStart)
      value$crlEnd <- as.numeric(value$crlEnd)
      
      if (length(value$bases) != length(value$qscores)) {
        # Truncate the longer vector to the shorter length so JS receives
        # matching-length arrays. This avoids rendering errors.
        min_len <- min(length(value$bases), length(value$qscores))
        warning(sprintf(
          "Length mismatch: bases (%d) vs qscores (%d) — truncating to %d",
          length(value$bases), length(value$qscores), min_len
        ))
        value$bases <- value$bases[seq_len(min_len)]
        value$qscores <- value$qscores[seq_len(min_len)]
      }

      value  # Return the validated and processed data
    }, error = function(e) {
      warning("Error processing sequence data: ", e$message)
      NULL
    })
  }
}