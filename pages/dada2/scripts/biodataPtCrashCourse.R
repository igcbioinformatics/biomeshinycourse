compareSampleNames <- function(fwdPathList, revPathList, splitSymbol, pickElement) {
  
  # `compareSamplesNames()` function takes as input two path lists
  # it retrieves the `basename()`, i.e., the last directory folder or file
  # it split it by `splitSymbol` and picks up `splitSymbol` one of the elements
  
  fwdNames <- sapply(X = strsplit(x = basename(path = fwdPathList), split = splitSymbol), FUN = `[`, pickElement)
  revNames <- sapply(X = strsplit(x = basename(path = revPathList), split = splitSymbol), FUN = `[`, pickElement)
  
  if (length(fwdNames) == length(revNames)) {
    for (name in 1:length(revNames)) {
      if (fwdNames[name] != revNames[name]) {
        #print("")
        print(paste0(fwdNames[name], " does not match ", revNames[name], "!"))
        print("`fwdPathList` and `revPathList` are not sorted!")
        #print("")
      }
      else {
        #print("")
        print(paste0(fwdPathList[name], " and ", revPathList[name], " are sorted!"))
        #print("")
      } 
    }
  }
  else {
    #print("")
    print("`fwdPathList` and `revPathList` have different length!")
    #print("")
    }
}



#-------------------------------------------------------------------------------------------------------------------------------------------

countUniqueFromDadaObjcList <- function(dadaObj) {
  
  ### count the total no. of unique seqs in the dada2 object

  denoise <- sapply(X = dadaObj, FUN = summary)[1,]
  denoiseOut <- as.numeric(denoise)
  names(denoiseOut) <- names(dadaObj)
  
  uniq <- sapply(X = dadaObj, FUN = summary)[11,]
  uniqOut <- as.numeric(uniq)
  names(uniqOut) <- names(dadaObj)
  
  listOut <- list()
  listOut[["unique"]] <- uniqOut
  listOut[["denoised"]] <- denoiseOut
  
  return(listOut)
  
}



#-------------------------------------------------------------------------------------------------------------------------------------------

countMergedSeqsFromDadaObjcList <- function(mergedObj) {
  
  ### count the total no. of merged seqs in the dada2 object
  
  merged <- sapply((sapply(X = mergedObj, `[`, 2)), sum)
  mergedOut <- as.numeric(merged)
  names(mergedOut) <- names(mergedObj)
  
  return(mergedOut)
  
}



#-------------------------------------------------------------------------------------------------------------------------------------------

convertTab2Biom <- function(inFile, outFile) {
  
  if (system("command -v biom", ignore.stdout = TRUE, ignore.stderr = TRUE) !=0)  {

    stop("biom program is not installed or it is not accessible!\n  Exiting...")
    
  }
  
  system(paste("biom convert", "-i", inFile, "-o", outFile, "--to-hdf5", 
               '--table-type="OTU table"', "--process-obs-metadata taxonomy"))
  
}



#-------------------------------------------------------------------------------------------------------------------------------------------

tax2biom <- function(taxTable) {
  
  lenNrow = nrow(taxTable)
  lenNcol = ncol(taxTable)
  id = c()
  fullTax = c()
  
  
  for (irow in 1:lenNrow) {
    rowID = as.character(rownames(taxTable)[irow])
    id = append(id, rowID)
    tax = ""
    
    for (icol in 1:lenNcol) {
      colName <- colnames(taxTable)[icol]
      colSign <- paste0(tolower(strtrim(colName,1)), "__")
      
      if (tax == "") {
        tax <- paste0(colSign, taxTable[irow,icol])
        #tax <- taxTable[irow,icol]
      }
      else {
        tax <- paste(tax, paste0(colSign, taxTable[irow,icol]), sep = "; ")
        #tax <- paste(tax, taxTable[irow,icol], sep = "; ")
      }
    }
    fullTax <- append(fullTax, tax)
  }
  
  taxTable2Biom <- data.frame(id, fullTax)
  colnames(taxTable2Biom) <- c("ASV_ID", "taxonomy")
  
  return(taxTable2Biom)
}




#-------------------------------------------------------------------------------------------------------------------------------------------
