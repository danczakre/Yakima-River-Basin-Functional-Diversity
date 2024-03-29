# Modified version of comdistnt leveraging Rfast
require(Rfast)

comdistnt_fast = function (comm, dis, abundance.weighted = FALSE, exclude.conspecifics = FALSE) 
{
  dat <- match.comm.dist(comm, dis)
  comm <- dat$comm
  dis <- dat$dist
  N <- dim(comm)[1]
  comm <- decostand(comm, method = "total", MARGIN = 1)
  comdisnt <- matrix(nrow = N, ncol = N)
  if(length(which(is.na(dis))) > 0){stop("Something is amiss with you cophenetic distances.")}
  for (i in 1:(N - 1)) {
    for (j in (i + 1):N) {
      sppInSample1 <- colnames(comm[i, comm[i, ] > 0, 
                                    drop = FALSE])
      sppInSample2 <- colnames(comm[j, comm[j, ] > 0, 
                                    drop = FALSE])
      if ((length(sppInSample1) >= 1) && (length(sppInSample2) >= 
                                          1)) {
        sample.dis <- dis[sppInSample1, sppInSample2, 
                          drop = FALSE]
        if (exclude.conspecifics) {
          sample.dis[sample.dis == 0] <- NA
        }
        
        sample1NT <- rowMins(sample.dis, value = T)
        sample1NT[sample1NT == Inf] <- NA
        names(sample1NT) = row.names(sample1NT)
        
        sample2NT <- colMins(sample.dis, value = T)
        sample2NT[sample2NT == Inf] <- NA
        names(sample2NT) = colnames(sample2NT)
        
        if (abundance.weighted) {
          sample1.weights <- as.numeric(comm[i, sppInSample1])
          sample2.weights <- as.numeric(comm[j, sppInSample2])
          if (any(is.na(sample1NT))) {
            miss <- which(is.na(sample1NT))
            sample1NT <- sample1NT[-miss]
            sample1.weights <- sample1.weights[-miss]
            sample1.weights <- sample1.weights/sum(sample1.weights)
          }
          if (any(is.na(sample2NT))) {
            miss <- which(is.na(sample2NT))
            sample2NT <- sample2NT[-miss]
            sample2.weights <- sample2.weights[-miss]
            sample2.weights <- sample2.weights/sum(sample2.weights)
          }
          sampleNT <- c(sample1NT, sample2NT)
          sample.weights <- c(sample1.weights, sample2.weights)
          comdisnt[i, j] <- weighted.mean(sampleNT, 
                                          sample.weights, na.rm = TRUE)
        }
        else {
          comdisnt[i, j] <- mean(c(sample1NT, sample2NT), 
                                 na.rm = TRUE)
        }
      }
      else {
        comdisnt[i, j] <- NA
      }
    }
  }
  rownames(comdisnt) <- colnames(comdisnt) <- rownames(comm)
  return(as.dist(t(comdisnt)))
}
