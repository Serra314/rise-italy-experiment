# RISE Italy experiment

Bayesian ETAS implementation using Inlabru for the next time-dependent models RISE forecasting experiment for Italy 

The project contains four folders:

1. **daily_forecast** : a folder where the produced daily forecast will be stored.

2. **input** : a folder where the data used to produce the forecast is stored. It contains the ISIDE catalogue (ISIDE_catalog_selected_PyCSEP), the testing region (testing.area.dat) and the collection region (collection.area.dat).

3. **source** : a folder where the code to produce the forecast is stored. It contains two files: source_inlabru_utilities.R containing the functions to fit the model and produce the forecast, and source_inlabru_1dayfore.R which is the script that needs to be runned to produce the forecast. 

4. **documentation** : a folder containing a document illustrating the implemented model (Theory_doc.pdf), a document giving more details on the approximation performed (Theory_supplement.pdf), a document listing all the R-packages needed to run the code (Code_doc.pdf) and containing instruction on how to run the code. For more details on the R-packages we have provided our complete session information as an R object (session.Info.Rds). 

The .txt file parameters contains imporant information for the script producing the forecast. It needs to contain the starting and end dates of the forecast, provided as forecastStartDate = %Y/%m/%d %h:%m:%s and forecastEndDate = %Y/%m/%d %h:%m:%s, and the path for the data used for model fitting, provided as PathCatalogData = ..path, and the path of the folder where the forecast will be stored, provided as PathOutputData = ..path (the two points before the actual path are mandatory). Please refer to the provided parameter document for an example.

To produce a forecast is sufficient to open a terminal in the present folder and run the following command

`R < source/source_inlabru_1dayfore.R --save`
