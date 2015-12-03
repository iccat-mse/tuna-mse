### Read in EffDIS data from ICCAT ####

library(RODBC)
# postgres server running locally (in gedt /etc/odbc.ini)
chan <- odbcConnect("mse-local", case="postgresql", believeNRows=FALSE)
sqlTables(chan)  #List all tables in the DB
mydata <- sqlFetch(chan, "some_table") # Return a table as a dataframe
odbcClose(chan)

#postgres server at tuna-cc1
chan <- odbcConnect("effdis-tuna-cc1", case="postgresql", believeNRows=FALSE)
sqlTables(chan)  #List all tables in the DB
mydata <- sqlFetch(chan, "t2ce") # Return a table as a dataframe
odbcClose(chan)
