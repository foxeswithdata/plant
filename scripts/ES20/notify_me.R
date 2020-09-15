
library("pushoverr")

send_to_elisa <- function(msg){
  pushover(message = msg, 
           user="uxcdpei2f3h4smz3wzd592anu8g5p1",
           app="ak3mwpq96ziwqm8kfsg69ogy7sb1y9")
}


