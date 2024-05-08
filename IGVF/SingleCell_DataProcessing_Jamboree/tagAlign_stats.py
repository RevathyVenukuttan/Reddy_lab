#!/usr/bin/env python
import os
import pandas as pd

def getManualAnnotationLabel(file):
    print("file is {}".format(file))
    Cell=file.split(subdir+"_")[-1].split(".tsv")[-2]
    print("Cell is {}".format(Cell))
    return(Cell)

def getCellClusterID(file,subdir):
    ClusterID="Cluster{}".format(index+1)
    print("ClusterID is {}".format(ClusterID))
    CellClusterID = subdir+"_"+ClusterID
    print("CellClusterID is {}".format(CellClusterID))
    return(CellClusterID)

def getMeanRNAUMIsPerCell():
    return()


clusters_folder_path = "/data/gersbachlab/Revathy/IGVF/Jamboree/SingleCellData/notebooks/clusters/"

# List all subdirectories under the "clusters" folder
subdirectories = [name for name in os.listdir(clusters_folder_path) if os.path.isdir(os.path.join(clusters_folder_path, name))]
# print(subdirectories)
# Create an empty dictionary to store the files for each subdirectory
subdir_files_dict = {}

# List all files under each subdirectory and store in the dictionary
# for subdir in subdirectories[0:1]:  # for debug
for subdir in subdirectories:  # all files
    subdir_path = os.path.join(clusters_folder_path, subdir)
    files = [file for file in os.listdir(subdir_path) if os.path.isfile(os.path.join(subdir_path, file))]
    subdir_files_dict[subdir] = files
# print(subdir_files_dict)
# Display the files for each subdirectory
for subdir, files in subdir_files_dict.items():
    print("Subdirectory: {}".format(subdir))
    stats_output_file = os.path.join(clusters_folder_path,subdir,'{}_ClusterMetadata.tsv'.format(subdir))
    print("stats_output_file is {}".format(stats_output_file))
    if os.path.exists(stats_output_file):
        print("stats_output_file {} exists".format(stats_output_file))
        continue
        
    print("Files:")
    print(files)
    files_in_subdir = len(files)
    print("files_in_subdir {} is {}".format(subdir,files_in_subdir))
    # Define the headers of the stats TSV file
    headers = ['tagAlignFile','CellClusterID','ManualAnnotationLabel','nCells','MeanRNAUMIsPerCell','MeanATACFragmentsPerCell']
    # Create an empty DataFrame with the headers
    df_stats = pd.DataFrame(columns=headers)
    # for index, file in enumerate(files[0:1]): # debug
    for index, file in enumerate(files): # all files
        # Read the TSV file into a pandas DataFrame
        df = pd.read_csv(os.path.join(clusters_folder_path,subdir,file), delimiter='\t', \
                         names = ['chr','start','stop','cell_id','reads','strand'])
        # Now you can work with the DataFrame 'df'
        print(df.head(2))  # Example: display the first few rows of the DataFrame
        print("number of rows (+/-) is {}".format(df.shape[0]))
        # Filter out rows with "-" in the last column
        df_positive_strand = df[~(df.iloc[:, -1] == '-')]
        # Now you can work with the filtered DataFrame 'df_filtered'
        print(df_positive_strand.head(2))  # Example: display the first few rows of the filtered DataFrame
        number_of_fragments = df_positive_strand.shape[0]
        print("number_of_fragments = number of rows (+) is {}".format(number_of_fragments))
        # CellClusterID	ManualAnnotationLabel	nCells	MeanRNAUMIsPerCell	MeanATACFragmentsPerCell
        # Xu2020_Cluster1	CD4+ T cell		152	3500			5000
        CellClusterID = getCellClusterID(file,subdir)
        print("CellClusterID is {}".format(CellClusterID))
        ManualAnnotationLabel = getManualAnnotationLabel(file)
        print("ManualAnnotationLabel is {}".format(ManualAnnotationLabel))
        nCells = df_positive_strand['cell_id'].nunique()
        print("nCells is {}".format(nCells))
        MeanRNAUMIsPerCell = getMeanRNAUMIsPerCell()
        MeanATACFragmentsPerCell = number_of_fragments // nCells
        print("MeanATACFragmentsPerCell is {}".format(MeanATACFragmentsPerCell))
        data_to_Add = [[CellClusterID,ManualAnnotationLabel,nCells,"",MeanATACFragmentsPerCell]]
        # Row to add as a dictionary
        new_row = {'tagAlignFile': file,
                   'CellClusterID': CellClusterID, 
                   'ManualAnnotationLabel': ManualAnnotationLabel,
                   'nCells':nCells, 
                   'MeanRNAUMIsPerCell':MeanRNAUMIsPerCell,
                   'MeanATACFragmentsPerCell':MeanATACFragmentsPerCell}
        # Add the new row using loc
        df_stats.loc[len(df_stats)] = new_row
        print(df_stats.head())
    # Save the DataFrame to a TSV file for each dataset
    df_stats.to_csv(stats_output_file, sep='\t', index=False)