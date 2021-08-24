#' Get a list of constraints given a taxonomic dataset
#' @description Internal used to generate tree constraints
#' @keywords internal
#' @export

getListConstraints <- function(dataset, targetColumns, byClades = F) {
  classification <- lapply(1:length(targetColumns) - 1, function(x) {
    q <- unique(unlist(dataset[targetColumns[x + 1]]))
    na <- lapply(q, function(y) {
      dataset[which(dataset[, targetColumns[x + 1]] == y)[1], targetColumns[x]]
    })
    names(q) <- na

    q2 <- lapply(unique(names(q)), function(z) {
      q <- q[names(q) == z]
      names(q) <- NULL
      q
    })
    names(q2) <- unique(names(q))
    q2
  })
  names(classification) <- targetColumns



  ## Combine all species per first level

  firstLev <- lapply(classification[length(classification)][[1]], function(x) {
    paste0("(", paste(unlist(x), collapse = ","), ")")
  })
  names(firstLev) <- names(classification[length(classification)][[1]])

  ## Combine all families per superfamily

  additionalLevels <- rev(names(classification))[-1]

  if (byClades == F) {
    for (x in 1:length(additionalLevels)) {
      nextLev <- if (x == 1) {
        firstLev
      } else {
        nextLev
      }
      targetGroups <-
        classification[names(classification) == additionalLevels[x]][[1]]
      if (x == 1) {
        nextLevNested <- list()
        for (i in 1:length(targetGroups)) {
          nextLevP <- nextLev[names(nextLev) == targetGroups[i]]
          nextLevNested[[i]] <- paste(nextLevP, collapse = ",")
        }
        names(nextLevNested) <- names(targetGroups)
        nextLev <- nextLevNested
        nextLev
      } else {
        nextLevNested <- list()
        for (i in 1:length(targetGroups)) {
          nextLevP <- nextLev[names(nextLev) == targetGroups[[i]]]
          nextLevNested[[i]] <- paste(nextLevP, collapse = ",")
        }
        names(nextLevNested) <- names(targetGroups)
        nextLev <- nextLevNested
        nextLev
      }
      nextLev
    }
  } else {
    nextLevgroups <- list()
    for (x in 1:length(additionalLevels)) {
      nextLev <- if (x == 1) {
        firstLev
      } else {
        nextLev
      }
      targetGroups <-
        classification[names(classification) == additionalLevels[x]][[1]]
      if (x == 1) {
        nextLevNested <- list()
        for (i in 1:length(targetGroups)) {
          nextLevP <- nextLev[names(nextLev) == targetGroups[i]]
          nextLevNested[[i]] <- paste(nextLevP, collapse = ",")
        }
        names(nextLevNested) <- names(targetGroups)
        nextLev <- nextLevNested
        nextLev
      } else {
        nextLevNested <- list()
        for (i in 1:length(targetGroups)) {
          nextLevP <- nextLev[names(nextLev) == targetGroups[[i]]]
          nextLevNested[[i]] <- paste(nextLevP, collapse = ",")
        }
        names(nextLevNested) <- names(targetGroups)
        nextLev <- nextLevNested
        nextLev
      }
      nextLevgroups[[x]] <- nextLev
    }
    names(nextLevgroups) <- additionalLevels
    unTa <- lapply(nextLevgroups, function(x) unlist(lapply(x, unlist)))
    names(unTa[[1]]) <- names(firstLev)
    groups <- unlist(unTa)
    names(groups) <- gsub("^.*\\.", "", names(groups))
    nextLev <- as.list(groups)
  }
  return(nextLev)
}
