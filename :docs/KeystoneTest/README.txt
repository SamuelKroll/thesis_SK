Quick Readme

1. In the directory "lists" you need the edge table of the network to check. Each edge will be checked for dependency on control organisms (partial correlations)
2. In "OPU_Tables" you need the OPU/OTU tables that were used to generate the network
3. In the "Mapfiles" you need the mapfiles
4. You can run the commands in "Make_Testtable_Bacteria.r" which will process the OTU tables and fill the folder "Output_Testtables" with tables of organism abundances (these are the control organisms for partial correlations)
5. Then you can run the commands in the script "Test_Edge_Dependencies_BacV5Endo.r". This will fill the folder "Output_EdgeDependencies" with the new networks based on the partial correlations.
6. The script "Dependency_Evaluation" will then summarize the data. 