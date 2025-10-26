#' Create a chromatogram output element
#' @param outputId The ID of the output element
#' @param width Width of the plot container (CSS units)
#' @param height Height of the plot container (CSS units)
#' @export
chromatogramOutput <- function(outputId, width = "100%", height = "250px") {
  htmltools::tagList(
    # Container div that will be found by JavaScript
    htmltools::div(
      id = outputId,
      class = "chromatogram-output",
      style = sprintf("width: %s; height: %s; overflow-x: auto;", width, height)
    )
  )
}

#' Render a chromatogram plot
#' @param expr Expression that returns a list with traces, bases, peakLocations, and qualityScores
#' @export
renderChromatogram <- function(expr) {
  # Create a function from the expression
  func <- shiny::exprToFunction(expr, env = parent.frame())
  
  shiny::createRenderFunction(
    func = func,
    function(value, session, name, ...) {
      if (is.null(value) || !is.list(value)) return(NULL)
      # Convert any matrices/arrays to lists for JSON serialization
      value$traces <- lapply(value$traces, as.numeric)
      value$bases <- as.character(value$bases)
      value$peakLocations <- as.numeric(value$peakLocations)
      value$qualityScores <- as.numeric(value$qualityScores)
      value
    }
  )
}