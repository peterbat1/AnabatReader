#
# AnabatTools: an R script to read, display and analyse binary files produced by
# the Anabat frequency-division bat recording system.
#    Copyright (C) 2013, 2017  Peter D. Wilson
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# See http://www.gnu.org/licenses/ for a copy of the license.
#
# Development of an Anabat data file reader was undertaken by
# Peter D. Wilson between 11 and 17 December 2013, and a major revision completed
# 10 January 2017.
#
# This R package implements and extends earlier software (now lost in unfortunate
# circumstances) written by me in Borland Delphi between 2003 and 2011.
#
# This code (like its lost predecessor) is based on the comprehensive notes on
# Anabat file formats and how to read them, including example C++ code, provided
# by Chris Corben on his website:
#      http://users.lmi.net/corben/fileform.htm#Anabat%20File%20Formats
# (last accessed 11 January 2017).
#

require(stringr)
require(bitops)

# Some constants:
maxSize <- 16384 # As a legacy of its origins on early PCs running MSDOS, Anabat files represent
# Anabat files represent 16kbyte chunks of data




#################################################
#' Internal function to extract text header information from a binary Anabat file
#'
#' @param theData A binary data file.
#'
#' @return A named list with elements holding extracted fields from the text header. Fields are: TAPE, DATE, LOC, SPECIES, SPEC, NOTE and NOTE1.
#' @export
#'
getTextHeader <- function(theData)
{
  # At some point in the past it appears that there was a bug in the Anabat
  # software so that null or zero bytes would be written into the text header
  # fields instead of ASCII 32 for a space. It only shows up in a few files for
  # example in the NSW Bat Call Library data set. This patch should fix this
  # problem. PDW 3 May 2014
  stuff <- as.raw(theData[7:281])
  stuff[which(stuff==0)] <- charToRaw(" ")

  TAPE <- str_trim(rawToChar(stuff[1:8]))
  DATE <- str_trim(rawToChar(stuff[9:16]))
  LOC <- str_trim(rawToChar(stuff[17:56]))
  SPECIES <- str_trim(rawToChar(stuff[57:106]))
  SPEC <- str_trim(rawToChar(stuff[107:122]))
  NOTE <- str_trim(rawToChar(stuff[123:196]))
  NOTE1 <- str_trim(rawToChar(stuff[197:275]))

  # Following lines are for parsing text header in the original data:
  #   TAPE <- str_trim(rawToChar(as.raw(zz[7:14])))
  #   DATE <- str_trim(rawToChar(as.raw(zz[15:22])))
  #   LOC <- str_trim(rawToChar(as.raw(zz[23:62])))
  #   SPECIES <- str_trim(rawToChar(as.raw(zz[63:112])))
  #   SPEC <- str_trim(rawToChar(as.raw(zz[113:128])))
  #   NOTE <- str_trim(rawToChar(as.raw(zz[129:202])))
  #   NOTE1 <- str_trim(rawToChar(as.raw(zz[203:281])))

  return(list(TAPE = TAPE, DATE = DATE, LOC = LOC, SPECIES = SPECIES,
              SPEC = SPEC, NOTE = NOTE, NOTE1 = NOTE1))
}


#################################################
#' Internal function to extract Date and Time from a Type 132 binary Anabat file
#'
#' @param fileType Integer indicating Anabat file type
#' @param p Pointer (integer) to start of identified byte string containing encoded date and time
#' @param theData A binary Anabat data file
#'
#' @return A named list with the following integer-valued elements: DAY, MONTH, YEAR, HOURS, MINS, SECS, HUNDS, and MICROS
#' @export
#'
getDateTime <- function(fileType, theData)
{
  if (fileType == 132)
  {
    # Get file date and time:
    p <- strtoi("0x120") + 1
    YEAR <- theData[p] + 256*theData[p + 1]
    MON <- theData[p + 2]
    DAY <- theData[p + 3]
    HOURS <- theData[p + 4]
    MINS <- theData[p + 5]
    SECS <- theData[p + 6]
    HUNDS <- theData[p + 7]
    MICROS <- theData[p + 8] + 256*theData[p + 9]
  }
  else
  {
    YEAR <- 0
    MON <- 0
    DAY <- 0
    HOURS <- 0
    MINS <- 0
    SECS <- 0
    HUNDS <- 0
    MICROS <- 0
  }

  return(list(DAY = DAY, MON = MON, YEAR = YEAR, HOURS = HOURS,
              MINS = MINS, SECS = SECS, HUNDS = HUNDS, MICROS = MICROS))
}


#################################################
#' Internal function to extract key paramaters from a binary Anabat file
#'
#' @param p Pointer (integer) to the start of the byte string containing paramaters to be extracted
#' @param theData A binary Anabat data file
#'
#' @return A named list with the following elements: RES1, DIVRAT, VRES, and timeFactor
#' @export
#'
getParams <- function(p, theData)
{
  RES1 <- theData[p + 2] + 256*theData[p + 3]

  if (RES1 != 25000)
  {
    timeFactor <- 25000/RES1
  }
  else
  {
    timeFactor <- 1
  }

  DIVRAT <- theData[p + 4]

  VRES <- theData[p + 5]

  return(list(RES1 = RES1, DIVRAT = DIVRAT, VRES = VRES, timeFactor = timeFactor))
}


#################################################
#' Internal function to extract time-frequency data from a binary Anabat file for Filetype 129
#'
#' @param p Pointer (integer) to the start of the block containing call data
#' @param params Parameter set as produced by function getParams
#' @param theData A binary Anabat data file. These data are further processed by function \code{\link{calcFreq}} to derive frequency data
#'
#' @return A named list: timeData, last.t, and showDot
#' @export
#'
getData129 <- function(p, params, theData)
{
  if ((params$RES1 > 60000) || (params$RES1 < 10000)) return(NULL)

  t <- 1
  s <- 0
  lastdif <- 0
  timeData <- rep(0, maxSize)
  showDot <- rep(0, maxSize)
  nBytes <- length(theData)
  time <- 0

  while (p < nBytes)
  {
    if (theData[p] < 128)
    {
      dif <- theData[p]
      lastdif <- lastdif + dif
      time <- time + floor(params$timeFactor*lastdif + 0.5)
      timeData[t] <- time
      t <- t + 1
      p <- p + 1
    }
    else
    {
      if (theData[p] > 248)
      {
        s <- theData[p] - 248
        showDot[t:t + s-1] <- 1
        p <- p + 1
      }
      else
      {
        nShift <- theData[p]
        nShift <- bitwShiftR(nShift, 3)
        nShift <- bitwAnd(nShift, as.integer("0x0f"))
        dif <- (bitwAnd(theData[p], as.integer("0x07")))*256 + theData[p + 1]
        if (nShift > 0) dif <- bitwShiftL(dif, nShift)
        lastdif <- dif
        time <- time + floor(params$timeFactor*lastdif + 0.5)
        timeData[t] <- time
        t <- t + 1
        p <- p + 2
      }
    }
  }

  return(list(timeData = timeData, last.t = t - 1, showDot = showDot))
}


#################################################
#' Internal function to extract time-frequency data from a binary Anabat file for Filetypes 130, 131 and 132
#'
#' @param p Pointer (integer) to the start of the block containing call data
#' @param params Parameter set as produced by function getParams
#' @param fileType Integer indicating the file type: 130, 131 or 132
#' @param theData A binary Anabat data file
#'
#' @return A named list: timeData, last.t, and showDot. These data are further processed by function \code{\link{calcFreq}} to derive frequency data
#' @export
#'
getData130 <- function(p, params, fileType, theData)
{
  time <- 0
  dif <- 0
  lastdif <- 0
  t <- 1
  s <- 0
  timeData <- rep(0, maxSize)
  showDot <- rep(2, maxSize)
  showDot[0] <- 0
  showDot[1] <- 1

  nBytes <- length(theData)

  if ((params$RES1 > 60000) || (params$RES1 < 10000)) return(NULL)

  while ((p < nBytes) && (t < maxSize))
  {
    if (theData[p] < 128)
    {
      dif <- theData[p]
      if (dif > 63) dif <- -1*(bitFlip(dif,bitWidth=6) + 1)
      lastdif <- lastdif + dif
      time <- time + floor(params$timeFactor*lastdif + 0.5)
      timeData[t] <- time
      t <- t + 1
      p <- p + 1
    }
    else
    {
      if (theData[p] >= 224 ) # Show status
      {
        if (fileType > 130)
        {
          # Filetpes 131 and 132
          if (p > nBytes) break # =
          c <- bitwAnd(theData[p], 3)
          s <- theData[p + 1]
          if ((t + s - 1) > 16384) s <- 16384 - t # limit index to arrays
          showDot[t:(t + s - 1)] <- c
          p <- p + 2
        }
        else
        {
          # Filetype 130
          s <- theData[p] - 224
          if ((t + s - 1) > 16383) s <- 16384 - t
          showDot[t:(t + s - 1)] <- c
          p <- p + 1
        }
      }
      else
      {
        if ((128 <= theData[p]) && (theData[p] <= 159))
        {
          if ((p+1) > nBytes) break # =
          dif <- 256*bitwAnd(theData[p], as.integer("0x1f")) + theData[p + 1]
          lastdif <- dif
          time <- time  + floor(params$timeFactor*lastdif + 0.5)
          timeData[t] <- time
          t <- t + 1
          p <- p + 2
        }
        else
        {
          if ((160 <= theData[p]) && (theData[p] <= 191))
          {
            if ((p+2) > nBytes) break # =
            dif  <- 256*256*bitwAnd(theData[p], as.integer("0x1f")) + 256*theData[p + 1] + theData[p + 2]
            lastdif <- dif
            time <- time + floor(params$timeFactor*lastdif + 0.5)
            timeData[t] <- time
            t <- t + 1
            p <- p + 3
            #break
          }
          else
          {
            if ((192 <= theData[p]) && (theData[p] <= 239))
            {
              if ((p+3) > nBytes) break # =
              dif <- 256*256*256*bitwAnd(theData[p], as.integer("0x1f")) + 256*256*theData[p + 1] + 256*theData[p + 2] + theData[p + 3]
              lastdif <- dif
              time <- time + floor(params$timeFactor*lastdif + 0.5)
              timeData[t] <- time
              t <- t + 1
              p <- p + 4
            }
          }
        }
      }
    }
  }

  return(list(timeData = timeData, last.t = t - 1, showDot = showDot))
}


#################################################
#' Internal function to compute frequency values
#'
#' @param params Parameter set as produced by function \code{\link{getParams}}
#' @param timeData Named list produced by either \code{\link{getData129}} or \code{\link{getData130}}
#' @param N Index (integer) to the last item in the timeData data set
#'
#' @return Named list: freq, showDot
#' @seealso \code{\link{getParams}}, \code{\link{getData129}}, \code{\link{getData130}}
#' @export
#'
calcFreq <- function(params, timeData, N)
{
  DIVRAT <- params$DIVRAT

  freq <- rep(0, length(timeData))
  showDot <- rep(0, length(timeData))
  t <- 3
  showDot[1] <- 0
  showDot[2] <- 1

  Tmin <- ceiling(DIVRAT*4)
  Tmax <- floor(DIVRAT*250)
  if (Tmin < 48) Tmin <- 48
  if (Tmax > 12589) Tmax <- 12589

  while(t <= N)
  {
    td <- timeData[t] - timeData[t - 2]
    if ((td >= Tmin) && (td <= Tmax))
    {
      freq[t] <- trunc(DIVRAT*1000000/td)
      showDot[t] <- 2
    }
    else
    {
      freq[t] <- 0
      showDot[t] <- 0
    }

    t <- t + 1
  }

  return(list(freq = freq, showDot = showDot))
}


#################################################
#' Compute statistics from data in an Anabat object
#'
#' @description Computes some basic statistics on the data supplied in parameter thisObj. Data can be restricuted to a rectangular window by specifying a start and end time, and independently, a lower frequency bound. Statistics are only computed on data remaining within the window. By default the entire data set is analysed.
#'
#' @param thisObj An object of class Anabat produced by the function \code{\link{readAnabat}}
#' @param timeSpan Character value indicating if statistics are to be computed on the "full" data sequence or a "partial" data sequence. Default = "full"
#' @param timeLimits If a partial sequence, then this 2-element numeric vector gives the start time (in seconds) and the end time (in seconds) used to "window" the sequence
#' @param cleanBelow A frequency filter value. All points less than or equal to the value (in Hertz) will be cleaned from the data before statistics are computed. A value of NULL (default) indicates no frequency filtering
#'
#' @return A named list with basic call statistics \emph{for the selected window}:
#' \tabular{ll}{
#'   meanCallDuration \tab Mean duration of a call (in seconds)\cr
#'   meanCallRate\tab Mean number of calls per second\cr
#'   meanCallGap\tab Mean time between calls (in seconds)\cr
#'   freqRange\tab Two-element numeric vector holding Minimum and Maximum frequencies (Hz)\cr
#'   meanFreq\tab Mean frequency (Hz)\cr
#'   medianFreq\tab Median frequency (Hz)\cr
#'   }
#' @export
#'
#' @examples
#' d <- readAnabat(system.file("extdata","Miau-ne.01#",package="AnabatTools"))
#'
#' # Basic stats for whole data set:
#' stats1 <- callStats(d)
#'
#' # Which is the same as the raw statistics contained in the object created by readAnabat:
#'
#' # Window on time axis:
#' stats2 <- callStats(d, startTime = 2.4, endTime = 4)
#'
#' # Window on the frequency axis using cleanBelow:
#' stats3 <- callStats(d, cleanBelow = 50000)
#'
#' # Combined window:
#' stats3 <- callStats(d, startTime = 2.4, endTime = 4, cleanBelow = 50000)


#################################################
callStats <- function(thisObj, timeSpan = "full",
                      timeLimits = NULL, cleanBelow = NULL)
{
  if (class(thisObj) == "Anabat")
  {
    if (timeSpan %in% c("full","partial"))
    {
      # Instantiate working data vectors
      fd <- thisObj$frequencyData
      td <- thisObj$timeData/1E06 # convert to seconds

      # Apply frequency filter
      if (!is.null(cleanBelow))
      {
        if (is.numeric(cleanBelow))
        {
          goodDataInd <- which(fd <= cleanBelow)
          if (length(goodDataInd > 0)) fd[goodDataInd] <- NA
        }
        else
        {
          message("ERROR: AnabatTools:callStats: Paramater cleanBelow must be numeric")
          return(NULL)
        }
      }

      if (timeSpan == "full")
      {
        goodPtInd <- 1:length(td)
      }
      else
      {
        if ((is.vector(timeLimits) && (length(timeLimits) == 2) && is.numeric(timeLimits)))
        {
          if (all(findInterval(timeLimits,range(td))))
          {
            goodPtInd <- intersect(which(td >= timeLimits[1]),
                                   which(td <= timeLimits[2]))

            if (length(goodPtInd) == 0)
            {
              message("ERROR: AnabatTools:callStats: Values passed in timeLimits leave no points available for stats")
              return(NULL)
            }
          }
          else
          {
            message("ERROR: AnabatTools:callStats: Paramater timeLimits are out of range for thisObj")
            return(NULL)
          }
        }
        else
        {
          message("ERROR: AnabatTools:callStats: Paramater timeLimits must be numeric vector of length 2")
          return(NULL)
        }
      }

      # OK to proceed with stats calcs now:
      freqNA <- which(is.na(fd[goodPtInd]))
      fd <- fd[goodPtInd]
      if (length(freqNA) > 0) fd[freqNA] <- 0
      fd[which(fd > 0)] <- 1

      td <- td[goodPtInd]

      zz <- rle(fd)

      onDots <- which(zz$values == 1)
      nCalls <- length(onDots)
      offDots <- which(zz$values == 0)

      endPt <- cumsum(zz$lengths)
      startPt <- endPt - zz$lengths + 1

      eventDuration <- td[endPt] - td[startPt]
      meanCallDuration <- mean(eventDuration[onDots])
      meanCallRate <- nCalls/(td[endPt[nCalls]] - td[startPt[1]])

      meanCallGap <- mean(eventDuration[offDots])

      validFreqVals <- thisObj$frequencyData[goodPtInd][which(fd == 1)]

      freqRange <- range(validFreqVals)
      meanFreq <- mean(validFreqVals)
      medianFreq <- median(validFreqVals)

      callStatsObj <- list(meanCallDuration = meanCallDuration,
                           meanCallRate = meanCallRate,
                           meanCallGap = meanCallGap,
                           freqRange = freqRange,
                           meanFreq = meanFreq,
                           medianFreq = medianFreq)
      class(callStatsObj) <- "AnabatCallStats"
      return(callStatsObj)
      # return(list(meanCallDuration = meanCallDuration,
      #             meanCallRate = meanCallRate,
      #             meanCallGap = meanCallGap,
      #             freqRange = freqRange,
      #             meanFreq = meanFreq,
      #             medianFreq = medianFreq))
    }
    else
    {
      message("ERROR: AnabatTools:callStats: Paramater timeSpan must be one of 'full' (= default) or 'partial'")
      return(NULL)
    }
  }
  else
  {
    message("ERROR: AnabatTools:callStats: Paramater thisObj must be class Anabat")
    return(NULL)
  }
}


#################################################
#' Read and parse an Anabat zero-crossing data file returning an S3 object containing the recovered data
#'
#' @description The workhorse function which reads and parses an Anabat zero-crossing data file and creates an S3 object to contain it. When a valid data file is loaded, the function also calls \code{\link{callStats}} to compute call statistics for the whole or raw data. These results are incorporated into the returned Anabat object.
#'
#' @param fileName Path to an Anabat data file
#'
#' @return An S3 object containing the data extracted from fileName
#' @export
#'
#' @examples
#' d <- readAnabat(system.file("extdata","Miau-ne.01#",package="AnabatTools"))
#'
#' # Some exploration and quality checking:
#' summary(d)
#' plot(d)
readAnabat <- function(fileName)
{
  AnabatObj <- list()
  class(AnabatObj) <- "Anabat"
  AnabatObj$filename <- fileName

  rawData <- readBin(fileName, what = "integer", n = 16384, size = 1, signed = F)

  AnabatObj$fileType <- rawData[4]

  # Set Pointer to parameter table:
  p <- rawData[1] + 256*rawData[2] + 1

  # Fetch text header:
  AnabatObj$textHeader <- getTextHeader(rawData)

  # Get other file  parameters:
  params <- getParams(p,rawData)
  AnabatObj$params <- params

  #fileDateTime <- getDateTime(AnabatObj$fileType,rawData)
  #AnabatObj$recordingDate

  # Set pointer to data block:
  p <- rawData[p] + 256*rawData[p + 1]

  if (AnabatObj$fileType == 129)
  {
    # Get data from old 129 file
    timeResult <- getData129(p, params, rawData)
  }
  else
  {
    # Get data from file types 130, 131 & 132
    timeResult <- getData130(p, params, AnabatObj$fileType, rawData)
  }

  freqResult <- calcFreq(params, timeResult$timeData, timeResult$last.t)
  freq <- freqResult$freq
  showDot <- freqResult$showDot

  badPts <- which(freq == 0)
  freq[badPts] <- NA

  AnabatObj$frequencyData <- freq[1:timeResult$last.t]
  AnabatObj$showDot <- showDot[1:timeResult$last.t]
  AnabatObj$timeData <- timeResult$timeData[1:timeResult$last.t]

  AnabatObj$rawCallStats <- callStats(AnabatObj)

  return(AnabatObj)
}


#################################################
#' Generic function to produce a summary of an Anabat object on the console
#'
#' @param AnabatData An object of class Anabat as produced by the function readAnabat
#'
#' @return NULL
#' @export
#'
#' @examples
#' d <- readAnabat("blah")
#' summary(d)
summary.Anabat <- function(AnabatData)
{
  if (class(AnabatData) == "Anabat")
  {
    cat("----------------------------------------------------------\n")
    cat("File:", AnabatData$filename, "\n")
    cat("File type:", AnabatData$fileType, "\n\n")
    cat("Header:\n")
    cat("    TAPE   : ", AnabatData$textHeader$TAPE, "\n")
    cat("    LOC    : ", AnabatData$textHeader$LOC, "\n")
    cat("    DATE   : ", AnabatData$textHeader$DATE, "\n")
    cat("    SPECIES: ", AnabatData$textHeader$SPECIES, "\n")
    cat("    SPEC   : ", AnabatData$textHeader$SPEC, "\n")
    cat("    NOTE   : ", AnabatData$textHeader$NOTE, "\n")
    cat("    NOTE1  : ", AnabatData$textHeader$NOTE1, "\n\n")
    cat("Paramters:\n")
    cat("    RES1      : ", AnabatData$params$RES1, "\n")
    cat("    DIVRAT    : ", AnabatData$params$DIVRAT, "\n")
    cat("    VRES      : ", AnabatData$params$VRES, "\n")
    cat("    timeFactor: ", AnabatData$params$timeFactor, "\n\n")
    cat("Number of data points (samples): ", length(AnabatData$frequencyData), "\n\n")
    cat("Raw call statistics:\n")
    cat("    Mean Call Duration: ",AnabatData$rawCallStats$meanCallDuration,"seconds\n")
    cat("    Mean Call Rate    : ",AnabatData$rawCallStats$meanCallRate,"calls per second\n")
    cat("    Mean Call Gap     : ",AnabatData$rawCallStats$meanCallGap,"seconds\n")
    cat("    Frequency Range   : ",AnabatData$rawCallStats$freqRange[1],"to",AnabatData$rawCallStats$freqRange[2],"Hz\n")
    cat("    Mean Frequency    : ",AnabatData$rawCallStats$meanFreq,"Hz\n")
    cat("    Median Frequency  : ",AnabatData$rawCallStats$medianFreq,"Hz\n")
    cat("----------------------------------------------------------\n")
  }
  else
  {
    message("ERROR: AnabatTools:summary.Anabat: Paramater AnabatData must be class Anabat")
  }
}


#################################################
#' Generic function to print summary data for an object of class Anabat
#'
#' @description This is convenience function which calls \code{\link{summary.Anabat}} and produces the same output to the console.
#'
#' @param AnabatData An object of class Anabat as produced by the function \code{\link{readAnabat}}
#'
#' @return NULL
#' @export
#'
#' @examples
#' d <- readAnabat("blah")
#' print(d)
print.Anabat <- function(AnabatData)
{
  if (class(AnabatData) == "Anabat")
  {
    summary(AnabatData)
  }
  else
  {
    message("ERROR: AnabatTools:print.Anabat: Paramater AnabatData must be class Anabat")
  }
}


#################################################
#' Generic plot function for Anabat objects
#'
#' @param AnabatData An object of class Anabat as produced by the function readAnabat
#' @param showLines Logical flag to indicate if horizontal grid lines should be produced on the plot. Default is TRUE.
#' @param showLineType Integer (default = 2, "dashed") indicating the line type to be used for horizontal grid lines. See help for par and option "lty"
#' @param lineColour Colour to be used for horizontal grid lines. This can be any valid value indexing a colour from the currently active palette.
#' @param maxFreq Maximum frequency to be used for the plot expressed in Hertz. That is, 35kHz would be entered as 30000. Default = 125000.
#' @param minFreq Minimum frequency to be used for the plot expressed in Hertz. That is, 5kHz would be entered as 5000. Default = 0.
#' @param cleanBelow Apply a simple frequency filter so that any points with frequency less than cleanBelow (in Hertz as described above) will be omitted from the plot.
#' @param startTime Time in seconds for the left margin of the plot. Default = 0 (ie start of the data sequence).
#' @param endTime Time in seconds for the right margin of the plot. Default = maximum time point in the data (ie end of the data sequence).
#' @param dotPch Numeric or character value specifying the plot symbol. See help for points function in the base graphics package for options. Default is ".".
#' @param dotCex Size of the specified plot symbol as a decimal fraction of the standard character size. Default = 0.3.
#' @param dotCol Colour to be used for the plot symbol. This can be any valid value (integer or character) indexing a colour from the currently active palette. Default = "blue" from the default R palette.
#'
#' @return NULL
#' @export
#'
#' @examples
#' d <- readAnabat(system.file("extdata","Miau-ne.01#",package="AnabatTools"))
#'
#' # Basic plot:
#' plot(d)
#'
#' # Change plot colour using a named colour from the default R palette:
#' plot(d, dotCol = "tomato")
#'
#' # Window on time axis:
#' plot(d, startTime = 2.4, endTime = 4)
#'
#' # Window on the frequency axis using cleanBelow:
#' plot(d, cleanBelow = 50000)
#'
#' # Combined window:
#' plot(d, startTime = 2.4, endTime = 4, cleanBelow = 50000)
plot.Anabat <- function(AnabatData,showLines=T,showLineType=2,lineColour="grey50",
                           maxFreq=125000,minFreq=0,cleanBelow=NULL,startTime=NULL,endTime=NULL,
                           dotPch=16,dotCex=0.3,dotCol="blue")
{
  if (class(AnabatData) == "Anabat")
  {
    plotFreqData <- AnabatData$frequencyData

    if (!is.null(cleanBelow))
    {
      if (is.numeric(cleanBelow))
      {
        goodDataInd <- which(plotFreqData <= cleanBelow)
        if (length(goodDataInd > 0)) plotFreqData[goodDataInd] <- NA
      }
      else
      {
        message("ERROR: AnabatTools:plot.Anabat: Paramater cleanBelow must be numeric")
        return(NULL)
      }
    }

    if (is.null(startTime)) startTime <- 0
    if (is.null(endTime)) endTime <- max(AnabatData$timeData/1E06)

    plot(AnabatData$timeData/1E06,plotFreqData,
         xlim=c(startTime,endTime),ylim=c(minFreq,maxFreq),
         pch=dotPch,cex=dotCex,col=dotCol,xlab="Time (seconds)",ylab="Frequency (Hz)") #,
         #main=AnabatData$textHeader$SPECIES)

    mtext(substitute(bolditalic(zz),list(zz = trimws(AnabatData$textHeader$SPECIES))),line=2)
    mtext(paste("Location:",AnabatData$textHeader$LOC,"  Date:",AnabatData$textHeader$DATE),cex = 0.9, col = "grey30",line = 1)

    if (showLines)
    {
      # Find out y-tick positions, but for now...
      nLines <- (maxFreq+10000) %/% 10000
      for (i in 1:nLines)
      {
        abline(h=(i-1)*10000,lty=showLineType,col=lineColour)
      }
    }
  }
  else
  {
    message("ERROR: AnabatTools:plot.Anabat: Paramater AnabatData must be class Anabat")
    return(NULL)
  }
}

