#________________________________________________________
# `.` is the magrittr placeholder in `%>%` pipelines, e.g. 
# `x %>% sub(" ", "", .)`. Declare it so R CMD check doesn't flag it as an 
# undefined global variable.
#________________________________________________________
utils::globalVariables(".")
