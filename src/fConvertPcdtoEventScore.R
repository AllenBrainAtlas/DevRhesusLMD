# Find equivalent developmental ages between species based on 
# Translating Time event scores (Workman et al. 2012)

ConvertPcdtoEventScore <- function(ageIn, 
                                   speciesIn="human", 
                                   speciesOut="eventScore") {
  # speciesOut can also be "eventScore"
  if(speciesIn=="eventScore"){ I1 = 0; I2 = 0;  }
  if(speciesIn=="human")     { I1 = 3.167; I2 = 3.72;  }
  if(speciesIn=="human_bc")     { I1 = 3.167; I2 = 3.72;  }
  if(speciesIn=="macaque")   { I1 = 3.27;  I2 = 2.413; }
  if(speciesIn=="rat")       { I1 = 2.31;  I2 = 1.705; }
  if(speciesIn=="mouse")     { I1 = 2.145; I2 = 1.894; }
  
  if(speciesIn=="eventScore") {
    eventScore <- ageIn
  } else {
    eventScore <- (log(2^ageIn)-I1)/I2  
  }
  
  if(speciesOut=="eventScore") return(eventScore)
  if(speciesOut=="human")    { R1 = 3.167; R2 = 3.72;  }
  if(speciesOut=="human_bc")    { R1 = 3.167; R2 = 3.72;  }
  if(speciesOut=="macaque")  { R1 = 3.27;  R2 = 2.413; }
  if(speciesOut=="rat")      { R1 = 2.31;  R2 = 1.705; }
  if(speciesOut=="mouse")    { R1 = 2.145; R2 = 1.894; }
  ageOut <- log2(exp(eventScore*R2+R1))
  return(ageOut)
}
