# --- Function to format bases as HTML with quality-based background colors ---
# This function is useful for displaying sequence data in data tables 
# (e.g., reactable/DT) where the background color indicates the quality score.
#
# Arguments:
# bases: A character vector of bases (e.g., c("A", "T", "G", "C")).
# qscores: A numeric vector of corresponding Phred quality scores, or a character 
#          vector of ASCII-encoded FASTQ quality scores (Phred+33).
# qscore_type: A string indicating the type of qscores provided. 
#              - "numeric_phred" (default): qscores is a vector of numeric Phred scores.
#              - "fastq_ascii": qscores is a vector of single-character ASCII codes.
#
# Returns: A single HTML string.
format_bases_as_html <- function(bases, qscores, qscore_type = "numeric") {
  
  if (length(bases) != length(qscores)) {
    # If the lengths don't match, return an error message
    return("Error: Base sequence and quality scores must have the same length.")
  }
  
  # --- 0. Quality Score Conversion ---
  # Convert FASTQ ASCII scores to numeric Phred scores if specified
  if (qscore_type == "fastq") {
    # Assuming Phred+33 encoding
    # Convert ASCII characters to their integer code, then subtract 33
    qscores_numeric <- as.integer(charToRaw(paste(qscores, collapse = ""))) - 33
    # Check if conversion resulted in correct length, if input was a character vector
    if (length(qscores_numeric) != length(qscores)) {
      stop("Error: FASTQ ASCII conversion failed. Check input format.")
    }
  } else {
    # Assume qscores are already numeric Phred scores
    qscores_numeric <- qscores
  }
  
  # 1. Define color mapping function based on Phred scores
  # Q >= 30: Green, Q >= 20: Yellow/Amber, Q < 20: Red/Pink
  get_color <- function(q) {
    if (q >= 20) {
      #return("#4CAF50") # High Quality (Green)
      return("rgba(76, 175, 80, 0.6);")
    } else if (q >= 15) {
      #return("#FFC107") # Medium Quality (Yellow/Amber)
      return("rgba(255, 193, 7, 0.6);")
    } else {
      #return("#F44336") # Low Quality (Red/Pink)
      return("rgba(244, 67, 54, 0.6);")
    }
  }
  
  # 2. Map Q-scores to background colors and create tooltip text
  
  background_colors <- sapply(qscores_numeric, get_color)
  
  # Create the data-tooltip attribute content (Position, Base, QV)
  base_positions <- seq_along(bases)
  # Round quality scores to integers for cleaner display in the tooltip
  qscores_rounded <- round(qscores_numeric, 0) 
  
  # Using data-tooltip attribute for custom CSS tooltip
  tooltip_content <- paste0(
    "Pos: ", base_positions, 
    " | Base: ", bases, 
    " | QV: ", qscores_rounded
  )
  
  # 3. Create a vector of <span> tags with inline CSS and DATA-TOOLTIP attribute
  # We add the class 'base-tooltip' for CSS targeting
  # Styles ensure monospaced font, black text, and compact spacing.
  html_spans <- paste0(
    '<span class="base-tooltip" data-tooltip="', tooltip_content, 
    '" style="background-color:', background_colors, 
    '; color: #000000; padding: 1px 0px; margin: 0; line-height: 1.5; font-family: monospace; font-weight: lighter; font-size: 1.0em;">', 
    bases, 
    '</span>'
  )
  
  # --- Formatting Logic for Blocks and Lines ---
  n <- length(html_spans)
  bases_html <- ""
  
  if (n > 0) {
    if (n > 1) {
      # 3a. Define separators vector (length n - 1)
      separators <- rep("", n - 1)
      
      # Indices 1:(n-1) correspond to the position *after* the base at that index.
      
      # 3b. Add spaces (&nbsp;) for blocks of 10 (e.g., after bases 10, 20, 30...)
      # Use &nbsp; for consistent, unbreakable block separation.
      indices_10 <- which((1:(n-1) %% 10 == 0) & (1:(n-1) %% 100 != 0))
      separators[indices_10] <- "&nbsp;"
      
      # 3c. Add <br> for lines of 100 (e.g., after bases 100, 200, 300...)
      indices_100 <- which(1:(n-1) %% 100 == 0)
      separators[indices_100] <- "<br>"
      
      # 3d. Interleave the bases and separators, then collapse to a single string
      # This creates a vector like [span1, sep1, span2, sep2, ..., spanN]
      bases_and_separators <- c(
        mapply(c, html_spans[-n], separators, SIMPLIFY = FALSE),
        list(html_spans[n])
      )
      bases_html <- paste(unlist(bases_and_separators), collapse = "")
      
    } else {
      # Single base case
      bases_html <- paste(html_spans, collapse = "")
    }
  }
  # --- End Formatting Logic ---
  
  # 4. Define CSS for the custom, non-delayed tooltip, now positioned to the RIGHT
  tooltip_css <- '
    /* Styles for the individual base span */
    .base-tooltip {
      position: relative; 
      cursor: default;
    }
    
    /* Tooltip text box - appears instantly on hover */
    .base-tooltip[data-tooltip]:hover::after {
      content: attr(data-tooltip);
      position: absolute;
      z-index: 10;
      top: 50%; /* Center vertically */
      left: 100%; /* Start position right of the span */
      /* Shift left by 50% of its width, then add 8px offset */
      transform: translateY(-50%) translateX(8px); 
      
      /* Styling */
      background-color: #333;
      color: #fff;
      padding: 4px 8px;
      border-radius: 4px;
      white-space: nowrap; /* Keep content on one line */
      font-size: 1em;
      opacity: 1; 
      pointer-events: none; 
      transition: opacity 0.1s; 
    }
    
    /* Tooltip arrow (optional, pointing left towards the base) */
    .base-tooltip[data-tooltip]:hover::before {
      content: "";
      position: absolute;
      top: 50%; /* Center vertically */
      left: 100%; /* Position at the right edge of the span */
      transform: translateY(-50%) translateX(3px); /* Move slightly into the text box */
      
      /* Arrow shape (Pointing Left) */
      border-width: 5px;
      border-style: solid;
      border-color: transparent #333 transparent transparent; 
      z-index: 11;
      opacity: 1;
      pointer-events: none;
      transition: opacity 0.1s;
    }
  '
  
  # 5. Wrap the formatted HTML content in a scrollable container including the new CSS
  # 'white-space: pre-wrap' ensures the <br> tags break the line correctly.
  full_html <- paste0(
    '<div style="white-space: pre-wrap; word-break: break-all; max-height: 350px; overflow-y: auto; border: 1px solid #ddd; padding: 8px; font-size: 0.9em;">',
    '<style>', tooltip_css, '</style>', # Inject CSS for custom tooltip
    bases_html,
    '</div>'
  )
  
  return(full_html)
}