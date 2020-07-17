library(RDCOMClient)
library(fs)

search.phrase <- 'Conjunction Reports'

save.fldr <- 'C:\\Users\\Harri\\OneDrive\\Documents\\centauri\\conjunction_analysis\\conj_data' # Set a root folder to save attachments into
print(save.fldr)

outlook_app <- COMCreate("Outlook.Application")
#Sys.sleep(40)
search <- outlook_app$AdvancedSearch(
  "Inbox",
  paste0("http://schemas.microsoft.com/mapi/proptag/0x0037001E ci_phrasematch '", search.phrase, "'")
)

Sys.sleep(40) # Wait some time to allow search to complete

results <- search[['Results']]

index = c(results[['Count']])


attachments.obj <- results[[index]][['attachments']]
  
 
      save.path <- file.path(save.fldr, attachments.obj[[1]][['FileName']]) # Set the save path
      print(save.path)
      attachments.obj[[1]]$SaveAsFile(save.path) # Save the attachment
  
  