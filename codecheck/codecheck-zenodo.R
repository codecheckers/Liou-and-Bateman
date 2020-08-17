## Zenodo deposit; see vignette("codecheck_overview.Rmd")
## https://github.com/codecheckers/notes/blob/master/codechecker-instructions.md
## Run this code from the codecheck/ folder.  But read it carefully
## rather than simply running the whole script.


library("codecheck")
library("zen4R")
library("yaml")

## Part one: create a new record ---------------------------------------
## You need do this only once for the certificate.

## To interact with the Zenodo API, you need to create a token.  This
## should not be shared, or stored in this script.  Uncomment one of
## the two following lines to get your token stored in `my_token`.
## my_token = system("pass show codechecker-token", intern=TRUE)
my_token = readLines("~/Documents/Zenodo-token.txt")


## make a connection to Zenodo API
zenodo <- ZenodoManager$new(token = my_token)


## This will generate a new record every time that you run it and 
## save the new record ID in the codecheck configuration file:
record = create_zenodo_record(zenodo)



## Part two: upload the metadata --------------------------------------
## You may run this several times if you edit the metadata.

metadata = read_yaml( "../codecheck.yml")

record = get_zenodo_record(metadata$report)
codecheck:::set_zenodo_metadata(zenodo, record, metadata)




## Part three: upload the certificate ---------------------------------
## If you have already uploaded the certificate once, you will need to
## delete it via the web page before uploading it again.

codecheck:::set_zenodo_certificate(zenodo, record, "codecheck.pdf")

## You may also create a ZIP archive of of any data or code files that
## you think should be included in the CODECHECK's record.



## Now go to zenodo and check the record (the URL is printed
## by set_zenodo_metadata() ) and then publish.
