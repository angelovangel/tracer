
library(RcppRoll)
library(seqinr)
library(dplyr)
library(stringr)
library(ggplot2)
library(ggnewscale)

# qscore_mean
qscore_mean <- function(qvs) {
  probs <- 10^(-qvs/10)
  meanprob <- mean(probs)
  -10*log10(meanprob)
}

# returns numeric
crl <- function(qscores, window_size, qval) {
  if (length(qscores) < window_size) {
    qscores_w <- rep(1, window_size)
  } else {
    qscores_w <- RcppRoll::roll_mean(qscores, n = window_size, na.rm = T)
  }
  rl <- rle(qscores_w >= qval)
  
  # if there are no rl20 then this is false and crl20 is set to 0, using just max() returns -Inf 
  if(any(rl$values)) {
    crl <- max(rl$lengths[rl$values], na.rm = TRUE)
    
    #start and end of max len T 
    
    trueind <- which(rl$values) # true indexes 
    maxtrueind <- trueind[which.max(rl$lengths[which(rl$values)])] # longest true index
    # stupid R 1-indexing
    crl_start <- ifelse(maxtrueind == 1, 1, cumsum(rl$lengths)[maxtrueind - 1])
    crl_end <- cumsum(rl$lengths)[maxtrueind]
    
  } else {
    crl <- 0
    crl_start <- 0
    crl_end <- 0
  }
  
  
  
  return(
    list(crl = crl, crl_start = crl_start, crl_end = crl_end)
  )
}

get_ab1 <- function(abfile) {
  obj <- sangerseqR::read.abif(filename = abfile)
   
  phredscores <- obj@data$PCON.1
  # use sangerseqR::read.abif only
  df <- tibble::tibble(
    #sample = obj@abifRawData@data$SMPL.1,
    sample = obj@data$SMPL.1,
    containerID = obj@data$CTID.1,
    #well = obj@abifRawData@data$TUBE.1,
    well = obj@data$TUBE.1,
    #rundate = paste0(obj@abifRawData@data$RUND.1$year, "-", obj@abifRawData@data$RUND.1$month, "-", obj@abifRawData@data$RUND.1$day),
    rundate = paste0(obj@data$RUND.1$year, "-", obj@data$RUND.1$month, "-", obj@data$RUND.1$day),
    #rawSeqLen = obj@QualityReport@rawSeqLength,
    seq = obj@data$PBAS.1 %>% str_split(''),
    qscores = list(phredscores),
    rawSeqLen = obj@data$PBAS.1 %>% str_split('') %>% unlist() %>% length(),
    
    #rawMeanQscore = obj@QualityReport@rawMeanQualityScore,
    rawMeanQscore = qscore_mean(obj@data$PCON.1),
    meanSNratio = mean(obj@data$`SN%.1`),
    
    #basesQ20 = sum(phredscores >= 20),
    basesQ20 = sum(obj@data$PCON.1 >= 20),
    #basesQ30 = sum(phredscores >= 30),
    #basesQ40 = sum(phredscores >= 40),
    
    crl20 = crl(phredscores, 20, 20)$crl,
    #polymerLot = obj@abifRawData@data$SMLt.1,
    polymerLot = obj@data$SMLt.1,
    polymer_expdate = obj@data$SMED.1,
    machine = obj@data$MCHN.1,
    instrument = obj@data$HCFG.3,
    capillary = obj@data$LANE.1,
    #len_to_detector = obj@abifRawData@data$LNTD.1,
    gel_type = obj@data$GTyp.1,
    analysis_prot = obj@data$APrN.1,
    data_coll_modfile = obj@data$MODF.1,
    dyeset_name = obj@data$DySN.1
    )
  
    # raw data
    df$data <- list(data = obj@data)
    return(
      df
    )
}

# --- Core Plotting Function ---
# This function encapsulates the logic for processing the ABIF file 
# and generating the ggplot object.
# NOTE: This version requires the 'ggnewscale' package.
# 
plot_abif_chromatogram <- function(rawdata, type = 'rawsignal') {
  
  
  if (length(rawdata) >= 0) {
    abif_data <- rawdata 
  } else {
    stop('rawdata not valid')
  }
  
  # DATA.1, DATA.2, DATA.3 and DATA.4 - raw signal
  # DATA.9, DATA.10, DATA.11 and DATA.12 - analysed signal
  if (type == 'rawsignal') {
    D1 <- abif_data$data$DATA.1
    D2 <- abif_data$data$DATA.2
    D3 <- abif_data$data$DATA.3
    D4 <- abif_data$data$DATA.4
  } else {
    D1 <- abif_data$data$DATA.9
    D2 <- abif_data$data$DATA.10
    D3 <- abif_data$data$DATA.11
    D4 <- abif_data$data$DATA.12
  }
  
  # Set a target max number of points for visualization
  MAX_POINTS_FOR_PLOT <- ifelse(type == 'rawsignal', 2000, 6000) 
  run_length <- length(D1)
  
  # Extract the raw signal trace data (A, C, G, T)
  # read FWO.1 to get which base to which DATA
  if (str_length(abif_data$data$FWO.1) == 4) {
    fwo <- str_split(abif_data$data$FWO.1, '') %>% unlist()
  } else {
    fwo <- c('G', 'A', 'T', 'C')
  }
  
  
  traces <- data.frame(
    time = 1:run_length,
    
    G = D1,  # A trace
    A = D2, # C trace (Note: DATA.10 here refers to the C trace, 
    # despite the generic name. DATA.11 is G, DATA.12 is T.)
    T = D3, # G trace
    C = D4  # T trace
  )
  colnames(traces) <- c('time', fwo)
  
  # --- Optimized Subsampling Logic ---
  if (run_length > MAX_POINTS_FOR_PLOT) {
    # Calculate step size to reduce the trace to approx MAX_POINTS_FOR_PLOT
    step_size <- ceiling(run_length / MAX_POINTS_FOR_PLOT)
    traces <- traces[seq(1, run_length, by = step_size), ]
  }
  # -----------------------------------
  
  # Convert data to long format for ggplot2
  traces_long <- reshape2::melt(traces, id.vars = "time", variable.name = "Base", value.name = "Intensity")
  
  
  # Extract base calls and their qv
  base_qv <- abif_data$data$PCON.1 # this is phred, not time
  base_locations <- abif_data$data$PLOC.1 # time + data$B1Pt.1 ?
  bases <- strsplit(abif_data$data$PBAS.1, "")[[1]]
  
  # FIX for differing vector lengths: Find the minimum length and truncate all vectors.
  min_len <- min(length(base_locations), length(bases), length(base_qv))
  
  # Check if base call data exists and align it
  if (length(base_locations) > 0 && length(bases) > 0 && length(base_qv) > 0) {
    base_calls <- data.frame(
      time = base_locations[1:min_len],
      base = bases[1:min_len],
      quality = base_qv[1:min_len] # Include the Quality Value
    )
    # Filter out base calls that are outside the trace length
    base_calls <- base_calls[base_calls$time <= run_length, ]
    
    
    # Scale the quality score (e.g., 0-60) to an alpha value (0-1).
    max_phred <- 60 
    
    base_calls$alpha_val <- pmin(base_calls$quality / max_phred, 1.0)
    # Set a minimum alpha floor (e.g., 0.3) so low-quality scores aren't invisible
    # base_calls$alpha_val <- pmax(base_calls$alpha_val, 0.3)
    
    # dd a new column for conditional bar coloring
    base_calls <- base_calls %>%
      dplyr::mutate(quality_color = dplyr::case_when(
        quality < 20   ~ "red",
        quality < 25  ~ "orange",
        TRUE          ~ "darkgrey"
      ))
    
  } else {
    # Create an empty dataframe if no base calls or QV are found
    base_calls <- data.frame(time = numeric(0), base = character(0), quality = numeric(0), alpha_val = numeric(0), quality_color = character(0))
  }
  
  # Define HIGH-VISIBILITY base colors (A=Green, C=Blue, G=Orange, T=Red)
  base_colors <- c("A" = "#00D100", "C" = "#0000FF", "G" = "black", "T" = "#FF0000") # Brighter versions
  
  # Create the base plot object
  p <- ggplot(traces_long, aes(x = time, y = Intensity, color = Base)) +
    #ggrastr::geom_point_rast(type = 'l', linewidth = 0.6, alpha = 0.5) +
    geom_line(linewidth = 0.6, alpha = 0.5) +
    scale_color_manual(values = base_colors) +
    theme_minimal(base_size = 14) + # Increased base_size for overall text
    theme(
      plot.title = element_blank(), 
      legend.position = "none",
      #panel.grid.minor = element_blank(), 
      panel.grid.major.y = element_blank(),
      plot.margin = margin(t = 5, r = 5, b = 0, l = 0, unit = "pt"), 
      axis.title.x = element_blank(),
      axis.title.y = element_text(margin = margin(t = 0, r = 5, b = 0, l = 0, unit = "pt")),
      axis.ticks.x = element_line(), # Ensure tick lines are visible
      axis.ticks.length = unit(0.3, "cm") # Increased tick length
    )
  
  # Add base call layers and modify axes IF base calls exist
  if ((nrow(base_calls) > 0) & (type != 'rawsignal')) {
    # Determine maximum intensity for positioning the labels
    max_intensity <- max(traces_long$Intensity)
    
    # Position for the base letter (A, C, G, T)
    base_label_y_pos <- max_intensity * 1.2
    # Define the y-range for the quality bars
    qv_bar_bottom_y_pos <- max_intensity * 1.0
    qv_bar_region_height <- max_intensity * 0.15
    
    # Add basecall index column for the new axis labels
    base_calls$index <- 1:nrow(base_calls)
    
    # Determine axis breaks and labels for the basecall index
    
    # --- Major Breaks (Labels) every 50 bases ---
    major_break_indices <- seq(0, nrow(base_calls), by = 50)
    axis_breaks <- base_calls$time[major_break_indices]
    axis_labels <- base_calls$index[major_break_indices]
    
    # --- Minor Breaks (Ticks) every 10 bases ---
    minor_break_indices <- seq(0, nrow(base_calls), by = 10)
    minor_breaks <- base_calls$time[minor_break_indices]
    # -------------------------------------------
    
    p <- p + 
      # Add base call labels (A, C, G, T)
      geom_text(
        data = base_calls,
        aes(x = time, y = base_label_y_pos, label = base), 
        color = base_colors[base_calls$base], 
        size = 4, #
        family = "mono", fontface = "bold",
        vjust = 0, 
        inherit.aes = FALSE 
      ) +
      
      # Start a new color scale for the quality bars
      ggnewscale::new_scale_color() +
      
      # Add Quality Value bars
      #ggrastr::geom_linerange_rast(
      geom_linerange(
        data = base_calls,
        aes(x = time, 
            ymin = qv_bar_bottom_y_pos, 
            ymax = qv_bar_bottom_y_pos + (quality / max_phred * qv_bar_region_height),
            #alpha = alpha_val,
            color = quality_color
        ), 
        linewidth = 1.3,
        inherit.aes = FALSE 
      ) +
      
      # Define the new scales
      scale_color_identity() +
      
      # Expand y-axis limits and set new x-axis labs and breaks
      coord_cartesian(ylim = c(0, max_intensity * 1.3)) +
      # labs(
      #   x = "Basecall Index",
      #   y = "Fluorescence Intensity",
      #   color = "Nucleotide"
      # ) +
      scale_x_continuous(
        breaks = axis_breaks,
        labels = axis_labels,
        minor_breaks = minor_breaks,
        expand = expansion(mult = c(0, 0.02))
      )
    
  } else {
    # Fallback for plots with no base calls: use original time axis
    p <- p + 
      labs(
        x = "Scan Number (Time)",
        y = "Fluorescence Intensity",
        color = "Nucleotide"
      ) +
      scale_x_continuous(expand = expansion(mult = c(0, 0.02)))
  }
  return(p)
}

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
format_bases_as_html <- function(bases, qscores, qscore_type = "numeric", crl_start, crl_end) {
  if ( any(!is.numeric(crl_start), !is.numeric(crl_end)) ) {
    return("Error: CRL start and end not numeric!")
  }
  
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
  # Q >= 20: Green, Q >= 15: Yellow/Amber, Q < 15: Red/Pink
  # transparency is passed as arg in order to be able to vary it
  get_color <- function(q, opacity = 0.6) {
    if (q >= 20) {
      return("#4CAF50") # High Quality (Green)
      #return(paste0("rgba(76, 175, 80,", opacity, ");"))
    } else if (q >= 10) {
      return("orange") # Medium Quality (Yellow/Amber)
      #return(paste0("rgba(255, 193, 7,", opacity, ");"))
    } else {
      return("#F44336") # Low Quality (Red/Pink)
      #return(paste0("rgba(244, 67, 54,", opacity, ");"))
    }
  }
  
  # 2. Map Q-scores to colors and create tooltip text
  # use seq_along to use index for crl
  
  colors <- sapply(seq_along(qscores_numeric), function(i){
    get_color(
      qscores_numeric[i], 
      opacity = ifelse(i > crl_start & i < crl_end, 1, 0.4)
    )
  }
  )
  
  # # Use grey background for seq outside crl
  bg_colors <- sapply(seq_along(qscores_numeric), function(i) {
    ifelse(i < crl_start | i > crl_end, "lightgrey", "")
  }
  )
  
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
    '" style="background-color:', bg_colors, 
    '; color:', colors,'; padding: 1px 0px; margin: 0; line-height: 1.5; font-family: monospace; font-weight: normal; font-size: 1.0em;">', 
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
    
    /* === Shadow the hovered base by changing its background === */
    .base-tooltip:hover {
      background-color: #dac586 !important; /* Light blue shadow */
    }
    /* ============================================================= */
    
    /* Tooltip text box - appears instantly on hover */
    .base-tooltip[data-tooltip]:hover::after {
      content: attr(data-tooltip);
      position: absolute;
      z-index: 10;
      top: 50%; /* Center vertically */
      left: 100%; /* Start position right of the span */
      /* Shift left by 50% of its width, then add 18px offset */
      transform: translateY(-50%) translateX(18px); 
      
      /* Styling */
      background-color: #333;
      color: #fff;
      padding: 4px 8px;
      border-radius: 4px;
      white-space: nowrap; /* Keep content on one line */
      font-size: 1em;
      opacity: 0.7; 
      pointer-events: none; 
      transition: opacity 0.1s; 
    }
    
    /* Tooltip arrow (optional, pointing left towards the base)
    .base-tooltip[data-tooltip]:hover::before {
      content: "";
      position: absolute;
      top: 50%; /* Center vertically */
      left: 100%; /* Position at the right edge of the span */
      transform: translateY(-50%) translateX(1px); /* Move slightly into the text box */
      
      /* Arrow shape (Pointing Left) */
      border-width: 5px;
      border-style: solid;
      border-color: transparent #333 transparent transparent; 
      z-index: 11;
      opacity: 1;
      pointer-events: none;
      transition: opacity 0.1s;
    } */
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

