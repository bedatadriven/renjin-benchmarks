count <- list(seq="ATG", count=0)

reads <- readRDS("reads.rds")
#reads <- sapply(1:1e6, function(x) paste(sample(c("A","T","C","G"), 36, replace = TRUE), collapse = ""))
#reads_1e7 <- sapply(1:1e7, function(x) paste(sample(c("A","T","C","G"), 36, replace = TRUE), collapse = ""))

seqCount <- function(e1, e2) {
  hasSeq <- grepl(e1[["seq"]], e2)
  if(hasSeq) e1[["count"]] <- e1[["count"]] + 1;
  e1
}

check_count <- function(object) {
  errors <- character()
  length_seq <- length(object@seq)
  if (length_seq != 1) {
    msg <- paste("Seq is length ", length_seq, ".  Should be 1", sep = "")
    errors <- c(errors, msg)
  }
  length_count <- length(object@count)
  if (length_count != 1) {
    msg <- paste("Name is length ", length_count, ".  Should be 1", sep = "")
    errors <- c(errors, msg)
  }
  if (!(object@count >= 0)) {
    msg <- paste("Count has value ", object@count, ". Should be > 0", sep = "")
    errors <- c(errors, msg)
  }
  if (length(errors) == 0) TRUE else errors
}


setClass("countS4", slots = list(seq = "character", count = "integer"))
setClass("countS4_valid", slots = list(seq = "character", count = "integer"), validity = check_count)
setClass("countS4_new", slots = list(seq = "character", count = "integer"))
setClass("countS4_new_valid", slots = list(seq = "character", count = "integer"), validity = check_count)


setGeneric("seqCountS4", function(count, read) standardGeneric("seqCountS4"))

setMethod("seqCountS4", representation(count = "countS4", read = "character"), function(count, read) {
  if(grepl(count@seq, read)) count@count <- count@count + 1L
  return(count)
})

setMethod("seqCountS4", representation(count = "countS4_valid", read = "character"), function(count, read) {
  if(grepl(count@seq, read)) count@count <- count@count + 1L
  return(count)
})

setMethod("seqCountS4", representation(count = "countS4_new", read = "character"), function(count, read) {
  new_count <- if(grepl(count@seq, read)) count@count + 1L else count@count
  return(new("countS4_new", seq = count@seq, count = new_count))
})

setMethod("seqCountS4", representation(count = "countS4_new_valid", read = "character"), function(count, read) {
  new_count <- if(grepl(count@seq, read)) count@count + 1L else count@count
  return(new("countS4_new_valid", seq = count@seq, count = new_count))
})

countS4a <- new("countS4", seq = "ATG", count = 0L)
countS4b <- new("countS4", seq = "ATG", count = 0L)
countS4aval <- new("countS4_valid", seq = "ATG", count = 0L)
countS4bval <- new("countS4_valid", seq = "ATG", count = 0L)
countS4anew <- new("countS4_new", seq = "ATG", count = 0L)
countS4bnew <- new("countS4_new", seq = "ATG", count = 0L)
countS4anewval <- new("countS4_new_valid", seq = "ATG", count = 0L)
countS4bnewval <- new("countS4_new_valid", seq = "ATG", count = 0L)


no_s4 <- system.time( for(i in 1:1e6) count <- seqCount(count, reads[i]), gcFirst = TRUE)

s4 <- system.time( for(i in 1:1e6) countS4a <- seqCountS4(countS4a, reads[i]), gcFirst = TRUE)

METHOD <- selectMethod("seqCountS4", c("countS4", "character"))
s4_preselect <- system.time( for(i in 1:1e6) countS4b <- METHOD(countS4b, reads[i]), gcFirst = TRUE)


s4_valid <- system.time( for(i in 1:1e6) countS4aval <- seqCountS4(countS4aval, reads[i]), gcFirst = TRUE)

METHOD <- selectMethod("seqCountS4", c("countS4_valid", "character"))
s4_valid_preselect <- system.time( for(i in 1:1e6) countS4bval <- METHOD(countS4bval, reads[i]), gcFirst = TRUE)


s4_new <- system.time( for(i in 1:1e6) countS4anew <- seqCountS4(countS4anew, reads[i]), gcFirst = TRUE)

METHOD <- selectMethod("seqCountS4", c("countS4_valid", "character"))
s4_new_preselect <- system.time( for(i in 1:1e6) countS4bnew <- seqCountS4(countS4bnew, reads[i]), gcFirst = TRUE)


s4_new_valid <- system.time( for(i in 1:1e6) countS4anewval <-     seqCountS4(countS4anewval, reads[i]) , gcFirst = TRUE)

METHOD <- selectMethod("seqCountS4", c("countS4_new_valid", "character"))
s4_new_valid_preselect <- system.time( for(i in 1:1e6) countS4bnewval <- METHOD(countS4bnewval, reads[i]), gcFirst = TRUE)


timings <- cbind(
  c(
	"Normal function call",
	"S4",
	"S4, preselection",
	"S4, validation",
	"S4, preselection, validation",
	"S4, new()",
	"S4, preselection, new()",
	"S4, validation, new()",
	"S4, preselection, validation, new()"),
  rbind(
	no_s4,
	s4,
	s4_preselect,
	s4_valid,
	s4_valid_preselect,
	s4_new,
	s4_new_preselect, 
	s4_new_valid, 
	s4_new_valid_preselect) 
  )

print(timings)


