# ============================================================================
# Functions shared by the figure scripts
#
# 1. Reading a raw event-log spreadsheet and turning it into visits. Used for
#    the 26 reference recordings (Figure 1) and for the RNAscope breeding
#    pairs (Figure 5), so that every event log in the paper is read with the
#    same rules. The rules are described at the top of
#    Figure1_ethogram_panels_11162026.R and in the reproducible report.
#
#      read_events(file, dir)     one row per In/Out event, typing slips flagged
#      pair_visits(events, file)  one row per visit: visitor, day, duration_min
#      read_visits(file, dir)     both steps, keeping visits of 0 to 24 hours
#
# 2. Scoring another cohort with the CarePC1 rule (Figures 3 and 5).
#
#      load_care_rule(rule_file)  returns care_pc1(), care_level() and the cutoff
# ============================================================================

library(readxl)

# ---- Read one spreadsheet into a clean event table --------------------------
# Each cell keeps its own type (col_types = "list"), so a column can be
# recognised by its contents rather than by its header or position.
share <- function(column, test) mean(vapply(column, function(cell) isTRUE(test(cell)), logical(1)))

is_status  <- function(x) is.character(x) && tolower(trimws(x)) %in% c("in", "out")
is_visitor <- function(x) is.character(x) && nchar(trimws(x)) <= 8 && grepl("^(m|fe|fa)", tolower(trimws(x)))
is_date    <- function(x) inherits(x, "POSIXct") && as.numeric(format(x, "%Y")) > 2000
is_time    <- function(x) inherits(x, "POSIXct") && as.numeric(format(x, "%Y")) < 1901   # Excel clock times
is_day     <- function(x) is.numeric(x) && !is.na(x) && x == round(x) && abs(x) < 100

find_columns <- function(sheet) {
  best <- function(test, skip = integer(0)) {
    scores <- vapply(sheet, share, numeric(1), test = test)
    scores[skip] <- -1
    unname(which.max(scores))
  }
  status <- best(is_status); visitor <- best(is_visitor, status)
  date   <- best(is_date);   time    <- best(is_time)
  day    <- best(is_day, c(status, visitor, date, time))
  c(day = day, date = date, time = time, visitor = visitor, status = status)
}

read_events <- function(file, dir) {
  sheet <- as.list(suppressMessages(read_excel(file.path(dir, file), sheet = 1,
                                               col_names = FALSE, col_types = "list")))
  col <- find_columns(sheet)
  keep <- vapply(sheet[[col[["status"]]]], is_status, logical(1)) &
          vapply(sheet[[col[["date"]]]], is_date, logical(1)) &
          vapply(sheet[[col[["time"]]]], is_time, logical(1))
  cell <- function(name) sheet[[col[[name]]]][keep]

  visitor <- vapply(cell("visitor"), function(x) {
    if (!is_visitor(x)) return(NA_character_)                  # blank or a note
    first <- tolower(substr(trimws(x), 1, 2))
    if (startsWith(first, "m")) "Male" else if (first == "fe") "Female" else "Family"
  }, character(1))

  date_part <- as.Date(do.call(c, cell("date")))
  year      <- format(date_part, "%Y")
  usual     <- names(which.max(table(year)))
  typo      <- year != usual
  if (any(typo)) {
    cat(sprintf("  %s: %d event(s) dated %s moved to %s\n", file, sum(typo),
                paste(unique(year[typo]), collapse = "/"), usual))
    date_part[typo] <- as.Date(paste0(usual, format(date_part[typo], "-%m-%d")))
  }
  clock   <- do.call(c, cell("time"))
  seconds <- as.numeric(format(clock, "%H", tz = "UTC")) * 3600 +
             as.numeric(format(clock, "%M", tz = "UTC")) * 60 +
             as.numeric(format(clock, "%S", tz = "UTC"))

  events <- data.frame(
    row     = which(keep),                                          # spreadsheet row number
    day     = vapply(cell("day"), function(x) if (is_day(x)) as.numeric(x) else NA_real_, numeric(1)),
    date    = as.numeric(date_part),
    when    = as.numeric(date_part) * 86400 + seconds,              # seconds on one continuous clock
    visitor = visitor,
    status  = tolower(trimws(unlist(cell("status")))),
    stringsAsFactors = FALSE)

  # fill blank Day cells from the date
  blank <- is.na(events$day)
  if (any(blank) && any(!blank)) {
    offset <- as.numeric(names(which.max(table(events$day - events$date))))
    events$day[blank] <- events$date[blank] + offset
    cat(sprintf("  %s: Day filled in from the date for %d event(s)\n", file, sum(blank)))
  }

  # flag rows logged more than 6 h earlier than a row above them with the same date
  events$out_of_sequence <- FALSE
  latest <- list()
  for (i in seq_len(nrow(events))) {
    key <- as.character(events$date[i])
    if (!is.null(latest[[key]]) && latest[[key]] - events$when[i] > 6 * 3600) {
      events$out_of_sequence[i] <- TRUE
    } else if (is.null(latest[[key]]) || events$when[i] > latest[[key]]) {
      latest[[key]] <- events$when[i]
    }
  }
  events
}

# ---- Pair each In with the next Out of the same visitor ---------------------
pair_visits <- function(events, file) {
  open   <- list()                             # the In event currently waiting, per visitor
  visits <- list()
  for (i in seq_len(nrow(events))) {
    who <- events$visitor[i]
    if (is.na(who)) {
      # a note in the visitor cell: an Out directly below an open In closes that visit
      above <- if (i > 1) events[i - 1, ] else NULL
      if (events$status[i] == "out" && !is.null(above) && above$status == "in" &&
          !is.na(above$visitor) && !is.null(open[[above$visitor]])) {
        who <- above$visitor
        cat(sprintf("  %s row %d: note in the visitor cell, read as %s Out\n", file, events$row[i], who))
      } else next
    }
    if (events$status[i] == "in") {
      open[[who]] <- events[i, ]
    } else if (!is.null(open[[who]])) {
      start    <- open[[who]]
      duration <- (events$when[i] - start$when) / 60
      if (start$out_of_sequence || events$out_of_sequence[i]) {
        # try each flagged row 12 h and 24 h later; keep the shortest positive duration
        shift_in  <- if (start$out_of_sequence)     c(720, 1440) else 0
        shift_out <- if (events$out_of_sequence[i]) c(720, 1440) else 0
        options   <- outer(shift_out, shift_in, function(o, s) duration + o - s)
        if (any(options > 0)) {
          fixed <- min(options[options > 0])
          cat(sprintf("  %s rows %d-%d (%s): %.0f min as typed, read as %.0f min\n",
                      file, start$row, events$row[i], who, duration, fixed))
          duration <- fixed
        }
      }
      visits[[length(visits) + 1]] <- data.frame(visitor = who, day = start$day, duration_min = duration)
      open[[who]] <- NULL
    }
  }
  do.call(rbind, visits)
}

# ---- Both steps: one spreadsheet -> its valid visits -------------------------
read_visits <- function(file, dir) {
  visits <- pair_visits(read_events(file, dir), file)
  visits[visits$duration_min >= 0 & visits$duration_min <= 24 * 60, ]
}

# ---- Scoring another cohort with the CarePC1 rule ---------------------------
# CarePC1 is defined once, on the 26 reference recordings. Figure2_11162026.R
# writes the rule to Supplementary/derived_tables/CarePC1_rule.csv: the mean and SD of coordinated
# visits per day and of coordinated minutes per day, the two loadings, and the
# High/Low cutoff. Every other cohort is scored with that same rule, so "High"
# and "Low" mean the same thing in every figure:
#
#   CarePC1 = loading1 * (coordinated visits per day  - mean1) / SD1
#           + loading2 * (coordinated minutes per day - mean2) / SD2
load_care_rule <- function(rule_file) {
  if (!file.exists(rule_file)) stop("Run Figure2_11162026.R first: it writes ", rule_file)
  rule   <- read.csv(rule_file, stringsAsFactors = FALSE)
  r_vis  <- rule[rule$Variable == "Visits_Family_perDay", ]
  r_min  <- rule[rule$Variable == "Minutes_Family_perDay", ]
  cutoff <- rule$Threshold[!is.na(rule$Threshold)][1]
  list(cutoff     = cutoff,
       care_pc1   = function(visits_per_day, minutes_per_day)
         r_vis$Loading * (visits_per_day  - r_vis$Mean) / r_vis$SD +
         r_min$Loading * (minutes_per_day - r_min$Mean) / r_min$SD,
       care_level = function(score) ifelse(score > cutoff, "High", "Low"))
}
