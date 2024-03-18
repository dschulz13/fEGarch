#'Shortcut to fitting an EGARCH-type Model to Observations
#'
#'@export
#'

`%><%` <- function(rt, spec) {
  egarchfit(rt = rt, spec = spec)
}

#'Pipe Operator
#'
#'@export
#'

`%>%` <- magrittr::`%>%`

