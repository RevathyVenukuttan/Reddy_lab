#!/usr/bin/env python
import os
import pandas as pd

full_fragments_and_cell_type_labels = [("/data/gersbachlab/Revathy/IGVF/Jamboree/SingleCellData/notebooks/fragment_files/syn52118183_syn52128237_GM12878_10XMultiome.atac.filter.fragments.hg38.tsv.gz/GM12878_10XMultiome.atac.filter.fragments.hg38.tsv.gz",
                                        "/data/gersbachlab/Revathy/IGVF/Jamboree/SingleCellData/results/GM12878_10xMultiome_cell_types_v1.tsv.gz")]

local_clusters_fld = os.path.join(os.getcwd(),"clusters")
os.makedirs(local_clusters_fld, exist_ok=True)
local_path_to_download = os.path.join(os.getcwd(),"fragment_files")

def remove_file(filename):
    print("remove_file method: {}".format(filename))
    if os.path.exists(filename):
        os.remove(filename)
        
import gzip
import pandas as pd

def read_tsv_gz_to_dataframe_skipping_comments_and_empty_lines(tsv_gz_file, comment_character='#'):
    # Open the compressed file using gzip
    print("tsv_gz_file is {}".format(tsv_gz_file))
    with gzip.open(tsv_gz_file, 'rt') as file:
        # Skip comment lines and empty lines, and load the remaining data into a DataFrame
        df = pd.read_csv(file, delimiter='\t', comment=comment_character, skip_blank_lines=True)

    return df

# # Specify the path to the TSV.gz file
# tsv_gz_file = '/path/to/file.tsv.gz'

# # Read the TSV.gz file and create a DataFrame
# dataframe = read_tsv_gz_to_dataframe(tsv_gz_file)

# # Display the DataFrame
# print(dataframe)

    
def split_fragment_line_string(string):
    # Remove newline characters
    string = string.replace("\n", "")

    # Splitting by tab character
    split_list = string.split("\t")

    # Splitting the word before the last one by underscore
    last_word = split_list[-2]
    split_word = last_word.split("_")

    # Inserting the split word before the last one in the list
    split_list.insert(-1, split_word[0])
    split_list.insert(-1, split_word[1])

    # Concatenating values at index 5 and index 4 with underscore
    concatenated_value = split_list[5] + "_" + split_list[4]
    split_list.append(concatenated_value)

    return split_list

def convert_fragment_line_to_tagAlign(r):
#     chr1	10007	10175	ENCSR023FME_GAAGGTTCAAAGTGTCAGTCAA	1
    rows_str = ""
    r_list = r.split("\t")
    # print("r_list is {}".format(r_list))
    row1 = []
    row2 = []

    row1.append(r_list[0])
    row1.append(r_list[1])
    row1.append(str(int(r_list[1]) + 1))
    row1.append(r_list[3])
    row1.append('1')
    row1.append('+')
    # print("convert_to_tagAlign: row 1 is: {}".format(row1))

    row2.append(r_list[0])
    row2.append(str(int(r_list[2]) - 1))
    row2.append(r_list[2])
    row2.append(r_list[3])
    row2.append('1')
    row2.append('-')
    # print("convert_to_tagAlign: row 2 is: {}".format(row2))
    
    rows_str = "\t".join(row1)+"\n"+"\t".join(row2)+"\n"
    # print("rows_str is {}".format(rows_str))
    return rows_str


for local_file_tuple in full_fragments_and_cell_type_labels:

    # Dictionary to store the output file handles with names
    output_handles = {}
    local_fragment_file = local_file_tuple[0]
    full_cell_types_annotation_file_path = local_file_tuple[1]
    print("!!!!!local_fragment_file is {}".format(local_fragment_file))
    file_atac_dataset_id = local_fragment_file.split("/")[10].split('.')[0]
    print("file_atac_dataset_id is {}".format(file_atac_dataset_id))
 
    # This full_cell_types_annotation_file_path was given as an input
    df_cell_types_for_atac_dataset = read_tsv_gz_to_dataframe_skipping_comments_and_empty_lines(full_cell_types_annotation_file_path)
    cell_type_id_names_for_atac_dataset = list(set(df_cell_types_for_atac_dataset['cell_type_id']))
    print("number of cell_type_id_names_for_atac_dataset is {}".format(len(cell_type_id_names_for_atac_dataset)))
    # List of output text files with corresponding names
    # [("output1.txt", "file_1"), ("output2.txt", "file_2"), ("output3.txt", "file_3")]
    output_tagAlign_files_with_names = [(os.path.join(local_clusters_fld,file_atac_dataset_id,"tagAlign_{}_{}.tsv".format(file_atac_dataset_id,cell_type_name_id)),
                                        cell_type_name_id) for cell_type_name_id in cell_type_id_names_for_atac_dataset]


    # this will make sure that we will not run the same tagAlign twice.
    tagAlign_exists = [os.path.exists(output_tagAlign_file[0]) for output_tagAlign_file in output_tagAlign_files_with_names]
    print("tagAlign_exists is {}".format(tagAlign_exists))
#     if sum(tagAlign_exists) >0:
#         print("output_tagAlign_files_with_names {} is at work or was already downloaded. continue".format(output_tagAlign_files_with_names))
#         continue # either started by annother process or already was processed
#     else:
    print("!!!output_tagAlign_files_with_names {}. open files".format(output_tagAlign_files_with_names))
    for tag_file_path, tag_file_cell_type_name in output_tagAlign_files_with_names:
#         print("tag_file_path is {}".format(tag_file_path))
#         print("os.path.dirname(tag_file_path) is {}".format(os.path.dirname(tag_file_path)))
        os.makedirs(os.path.dirname(tag_file_path), exist_ok=True)
        output_handles[tag_file_cell_type_name] = open(tag_file_path, "w")

    # print("!!!output_tagAlign_files_with_names {}".format(output_tagAlign_files_with_names))
    print("open local_fragment_file {}".format(local_fragment_file))
    with gzip.open(local_fragment_file, "rt") as infile:
        missing_bc = 0
        # Open the output files and store their handles in the list
        num_of_lines_written=0
        for line_number, line in enumerate(infile, start=1):
            # debug
#             if line_number > 20:
#                 continue 

            # here cases where bc_datasetId or datasetId_bc are being mixed between the fragments and the
            # cell type are being address. you can select the righ out_list for your experiment
            out_list = split_fragment_line_string(line)
            bc = out_list[-1]


#             Austin output: chrom, start, end, bc, rem = line.rstrip('\n').split('\t', 5)
            out_line_to_print = "{}\t{}\t{}\t{}\t{}\n".format(out_list[0],out_list[1],out_list[2],bc,out_list[-2])
#             chr1	10007	10175	ENCSR023FME#ENCSR023FME_GAAGGTTCAAAGTGTCAGTCAA	1
            num_of_lines_written +=1
            returnTagAlign = convert_fragment_line_to_tagAlign(out_line_to_print)

            # write to the relevant cell type file
            # print("df_cell_types_for_atac_dataset[cell_id][0:5] is {}".format(df_cell_types_for_atac_dataset['cell_id'][0:5]))
            bc_exists_in_cell_type_atac_dataset = df_cell_types_for_atac_dataset[df_cell_types_for_atac_dataset['cell_id']==bc]
            # print("len(bc_exists_in_cell_type_atac_dataset) is {} for bc {}".format(len(bc_exists_in_cell_type_atac_dataset),bc))
            if len(bc_exists_in_cell_type_atac_dataset) ==1:
                tag_file_cell_type_id = df_cell_types_for_atac_dataset.loc[df_cell_types_for_atac_dataset['cell_id'] == bc, 'cell_type_id'].iloc[0]
                # print("tag_file_cell_type_id is {}".format(tag_file_cell_type_id))
                output_handles[tag_file_cell_type_id].write(returnTagAlign)
            else:
                missing_bc+=1

        for tag_file_path, tag_file_cell_type_name in output_tagAlign_files_with_names:
            print("tag_file_path is {}".format(tag_file_path))
            output_handles[tag_file_cell_type_name].close()    
        print("finished clustering local_fragment_file {} by cell type. for types {}".format(local_fragment_file, cell_type_id_names_for_atac_dataset))
        print("total missing bc are {}".format(missing_bc))