# read ab1 file
# output is a data table with quality stats

# usage
# purrr::map_dfr(file_list, get_ab1)

# or for large number of files use multi cores
# plan(multisession, workers = 9)
# furrr::future_map_dfr(august, get_ab1_quality)
library(RcppRoll)
library(sangeranalyseR)
library(dplyr)

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
  
  obj <- sangeranalyseR::SangerRead(readFileName = abfile, readFeature = 'Forward Read')
    if(obj@objectResults@creationResult) {
      
      # calculate bases >Q20, >Q30, >Q40 as they are not explicitly stored in the SangerRead S4 instance
      phredscores <- obj@QualityReport@qualityPhredScores
      
      # continuous read length (CRL) - The longest uninterrupted stretch of bases with a running QV average of 20 or higher
      # use rle
      # use window of 20 bp for Q scores to adhere to definition used in Sanger Analysis software
      if (length(obj@QualityReport@qualityPhredScores) < 20) {
        qscores_w20 <- rep(1, 20)
      } else {
        qscores_w20 <- RcppRoll::roll_mean(obj@QualityReport@qualityPhredScores, n = 20, na.rm = T)
      }
      #qscores_w20 <- RcppRoll::roll_mean(obj@QualityReport@qualityPhredScores, n = 20, na.rm = T)
      rl20 <- rle(qscores_w20 >= 20)
      #rl30 <- rle(qscores_w20 >= 30)
      # if there are no rl20 then this is false and crl20 is set to 0, using just max() returns -Inf 
      if(any(rl20$values)) {
        crl20 <- max(rl20$lengths[rl20$values], na.rm = TRUE)
      } else {
        crl20 <- 0
      }
      
      #crl30 <- max(rl30$lengths[rl30$values], na.rm = TRUE)
      
      #============= Raw signal =======================================
      # raw signal intensities, to be able to present them in a table as sparklines:
      # https://glin.github.io/reactable/articles/examples.html#embedding-html-widgets
      # store in the dataframe as list of values
      # the values are roll means to reduce number of points to around 100
      
      # signal_counts <- obj@abifRawData@directory@numelements[match("DATA.1", names(obj@abifRawData@data))]
      # signal <- rowMeans(
      #   cbind(
      #     RcppRoll::roll_max(obj@abifRawData@data$DATA.1, n = 100, by = 500), 
      #     RcppRoll::roll_max(obj@abifRawData@data$DATA.2, n = 100, by = 500), 
      #     RcppRoll::roll_max(obj@abifRawData@data$DATA.3, n = 100, by = 500), 
      #     RcppRoll::roll_max(obj@abifRawData@data$DATA.4, n = 100, by = 500)
      #     )
      #   ) %>%
      #   #RcppRoll::roll_max(n = 500, by = 500) %>%
      #   round(digits = 0)
      # #============= Raw signal =======================================
      
      df <- tibble::tibble(
        sample = obj@abifRawData@data$SMPL.1,
        well = obj@abifRawData@data$TUBE.1,
        rundate = paste0(obj@abifRawData@data$RUND.1$year, "-", obj@abifRawData@data$RUND.1$month, "-", obj@abifRawData@data$RUND.1$day),
        rawSeqLen = obj@QualityReport@rawSeqLength, 
        trimSeqLen = obj@QualityReport@trimmedSeqLength,
        trimStart = obj@QualityReport@trimmedStartPos,
        trimEnd = obj@QualityReport@trimmedFinishPos,
        rawMeanQscore = obj@QualityReport@rawMeanQualityScore,
        trimMeanQscore = obj@QualityReport@trimmedMeanQualityScore,
        remainingRatio = obj@QualityReport@remainingRatio,
        basesQ20 = sum(phredscores >= 20),
        basesQ30 = sum(phredscores >= 30),
        basesQ40 = sum(phredscores >= 40),
        crl20 = crl20,
        polymerLot = obj@abifRawData@data$SMLt.1,
        polymer_expdate = obj@abifRawData@data$SMED.1,
        machine = obj@abifRawData@data$MCHN.1,
        instrument = obj@abifRawData@data$HCFG.3,
        capillary = obj@abifRawData@data$LANE.1,
        len_to_detector = obj@abifRawData@data$LNTD.1,
        gel_type = obj@abifRawData@data$GTyp.1,
        analysis_prot = obj@abifRawData@data$APrN.1,
        data_coll_modfile = obj@abifRawData@data$MODF.1,
        dyeset_name = obj@abifRawData@data$DySN.1,
        
        # signal = list(signal) does not work as expected
        #crl30 = crl30
        )
      # df$signal <- list(signal) # add signal as list to df 
      # so you can do reactable(columns = list(signal = colDef(cell = function(values) {sparkline(str_split(values, "\\|") %>% unlist() %>% as.numeric(), type = 'line')} )))
      df$seq <- as.character(obj@primarySeq) %>% str_split("")
      df$seqtrimmed <- str_sub(as.character(obj@primarySeq), obj@QualityReport@trimmedStartPos, obj@QualityReport@trimmedFinishPos) %>% str_split("")
      df$qscores <- list(phredscores)
      # raw data
      df$data <- list(data = obj@abifRawData@data)
      return(
        df
      )
    } else {
      stop(obj@objectResults@errorMessage)
    }
}

# --- Core Plotting Function ---
# This function encapsulates the logic for processing the ABIF file 
# and generating the ggplot object.
# rawdata is the 
plot_abif_chromatogram <- function(rawdata) {
  
  
  if (length(rawdata) >= 0) {
    abif_data <- rawdata 
  } else {
    stop('rawdata not valid')
  }
  
  
  # Extract the raw signal trace data (A, C, G, T)
  run_length <- length(abif_data$data$DATA.9)
  
  traces <- data.frame(
    time = 1:run_length,
    # read FWO.1 to get which base to which DATA
    G = abif_data$data$DATA.9,  # A trace
    A = abif_data$data$DATA.10, # C trace (Note: DATA.10 here refers to the C trace, 
    # despite the generic name. DATA.11 is G, DATA.12 is T.)
    T = abif_data$data$DATA.11, # G trace
    C = abif_data$data$DATA.12  # T trace
  )
  
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
    
    # 🌟 QV to Alpha Mapping Logic 🌟
    # Scale the quality score (e.g., 0-60) to an alpha value (0-1).
    max_phred <- 60 
    
    base_calls$alpha_val <- pmin(base_calls$quality / max_phred, 1.0)
    # Set a minimum alpha floor (e.g., 0.3) so low-quality scores aren't invisible
    base_calls$alpha_val <- pmax(base_calls$alpha_val, 0.3) 
    
  } else {
    # Create an empty dataframe if no base calls or QV are found
    base_calls <- data.frame(time = numeric(0), base = character(0), quality = numeric(0), alpha_val = numeric(0))
  }
  
  # Define colors for the traces
  base_colors <- c("A" = "green", "C" = "blue", "G" = "black", "T" = "red")
  
  # Create the plot
  p <- ggplot(traces_long, aes(x = time, y = Intensity, color = Base)) +
    geom_line(linewidth = 0.6, alpha = 0.5) +
    scale_color_manual(values = base_colors) +
    labs(
      title = "",
      x = "Scan Number (Time)",
      y = "Fluorescence Intensity",
      color = "Nucleotide"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_blank(),
      legend.position = "none",
      panel.grid.minor = element_blank(), 
      # *** START FIX FOR Y-AXIS SPACE ***
      plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"), 
      axis.title.y = element_text(margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")),
      axis.ticks.length.y = unit(0, "pt"), 
      axis.title.y.left = element_text(margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")),
      axis.ticks.length.x = unit(0, "pt")
      # *** END FIX FOR Y-AXIS SPACE ***
    )
  
  # Add base call labels (with alpha) and Quality Value scores (as solid text)
  if (nrow(base_calls) > 0) {
    # Determine maximum intensity for positioning the labels
    max_intensity <- max(traces_long$Intensity)
    
    # Position for the base letter (A, C, G, T) - higher
    base_label_y_pos <- max_intensity * 1
    
    # Position for the Quality Value score (e.g., 40, 55, 12) - lower
    qv_label_y_pos <- max_intensity * 0.95 
    
    # 1. Add base call labels (A, C, G, T) with ALPHA MAPPING
    p <- p + geom_text(
      data = base_calls,
      aes(x = time, y = base_label_y_pos, label = base), #alpha = alpha_val), 
      # Use the base color for the text labels
      color = base_colors[base_calls$base], 
      size = 3.5,
      family = "mono",
      vjust = 0, 
      inherit.aes = FALSE 
    ) +
      
      # 2. Add Quality Value scores (e.g., 40) - SOLID BLACK TEXT
      geom_text(
        data = base_calls,
        # Convert quality score to character for the label
        aes(x = time, y = qv_label_y_pos, label = as.character(quality), alpha = alpha_val), 
        color = "black", # Solid black color
        size = 1.3,      # Smaller font size
        family = "sans", 
        vjust = 0, 
        inherit.aes = FALSE 
      ) +
      
      # Crucial: Interpret alpha_val directly as the alpha level
      scale_alpha_identity() + 
      # Expand y-axis limits to ensure both sets of labels are visible
      coord_cartesian(ylim = c(0, max_intensity * 1.25)) +
      scale_x_continuous(expand = expansion(mult = c(0, 0.02)))
    
  } else {
    p <- p + scale_x_continuous(expand = expansion(mult = c(0, 0.02)))
  }
  
  return(p)
}
