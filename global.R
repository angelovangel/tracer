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
      
      signal_counts <- obj@abifRawData@directory@numelements[match("DATA.1", names(obj@abifRawData@data))]
      signal <- rowMeans(
        cbind(
          RcppRoll::roll_max(obj@abifRawData@data$DATA.1, n = 100, by = 500), 
          RcppRoll::roll_max(obj@abifRawData@data$DATA.2, n = 100, by = 500), 
          RcppRoll::roll_max(obj@abifRawData@data$DATA.3, n = 100, by = 500), 
          RcppRoll::roll_max(obj@abifRawData@data$DATA.4, n = 100, by = 500)
          )
        ) %>%
        #RcppRoll::roll_max(n = 500, by = 500) %>%
        round(digits = 0)
      #============= Raw signal =======================================
      
      df <- tibble::tibble(
        sample = obj@abifRawData@data$SMPL.1,
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
        machine = obj@abifRawData@data$MCHN.1
        # signal = list(signal) does not work as expected
        #crl30 = crl30
        )
      df$signal <- list(signal) # add signal as list to df 
      # so you can do reactable(columns = list(signal = colDef(cell = function(values) {sparkline(str_split(values, "\\|") %>% unlist() %>% as.numeric(), type = 'line')} )))
      df$seq <- as.character(obj@primarySeq) %>% str_split("")
      df$seqtrimmed <- str_sub(as.character(obj@primarySeq), obj@QualityReport@trimmedStartPos, obj@QualityReport@trimmedFinishPos) %>% str_split("")
      df$qscores <- list(phredscores)
      df
      
    } else {
      stop(obj@objectResults@errorMessage)
    }
}
