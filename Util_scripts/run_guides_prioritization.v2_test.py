# #!/usr/bin/env python
import numpy as np
import pandas as pd
import argparse

def overlapping(a, b, maximum_overlap=0):
    try:
        return not ((a['target site start coordinate']+maximum_overlap > b['target site end coordinate']) or \
                   (a['target site end coordinate']-maximum_overlap < b['target site start coordinate']))
    except:
        return not ((a['target site start coordinate'][0]+maximum_overlap > b['target site end coordinate']) or \
                   (a['target site end coordinate'][0]-maximum_overlap < b['target site start coordinate']))

def sort_grnas(x, maximum_overlap=5):
    selected = []
    for ix, i in x.iterrows():
        if not np.any([overlapping(x.loc[ii, :], i, maximum_overlap=maximum_overlap) for ii in selected]):
            selected.append(ix)
    return x.loc[selected, :]

def main():
    parser = argparse.ArgumentParser("""
    Prioritize gRNAs from GuideScan given their accessibility, while controling for the 
    maximum overlap among any two guides in the pool of selected guide RNAs.
    """)
    parser.add_argument('-guidescan-file', required = True, type=open, 
                        help="""Guidescan v1.0 batch results""")
    parser.add_argument('-outfile', required = True, type=str, 
                        help="""Guidescan results with overlap-controlled guides
                        ranked per region (in decreasing order of importance)""")
    parser.add_argument('--guides-accesibility-coverage', required=False,
                       help = """Accesibility counts table as computed by deeptools bamCoverage.
                       The first 3 columns are in BED3 format, followed by count columns. If more
                       than one count column, the mean will be taken. If the count columns should
                       be grouped, use the --library-sizes argument, one per group, to first compute
                       CPMs and then to take group means.""")
    parser.add_argument('--library-sizes', type=str, nargs='+', action='append', 
                        help="""List of library sizes in the order of the count columns in the 
                        coverage file. It is used to compute counts per millions.
                        The values can be provided directly or in files.""",
                        required=False)
    parser.add_argument('--max-overlap', type=int, default=5,
                        help="""Maximum overlap allowed between gRNAs selected.""")

    args = parser.parse_args()

    # Load GuideScan gRNA data
    headers = ['chromosome','target site start coordinate','target site end coordinate',
               'gRNA','cutting efficiency score','cutting specificity score','strand',
               'offtargets sum','offtargets summary','annotation','gRNA label']
    guide_scan_df = pd.read_csv(args.guidescan_file, names = headers, sep=',')
    guide_scan_df.index = guide_scan_df.apply(lambda x: "%s_%d_%d" % (x[0], 
                                                                      x[1],#-1 if x['strand'] == '+' else x[1]+2, 
                                                                      x[2]),#-2 if x['strand'] == '+' else x[2]+1), 
                                              axis=1)

    # If accessibility data provided, rank guides by it
    if args.guides_accesibility_coverage:
        counts_tmp = pd.read_csv(args.guides_accesibility_coverage, sep="\t", 
                                 comment="#", header=None)
        counts_tmp.index = counts_tmp.apply(lambda x: "%s_%d_%d" % \
                                            (x[0], x[1], x[2]), axis=1)
        counts_tmp = counts_tmp.drop_duplicates()
        print(counts_tmp.head())
    
        # Load library sizes and take the mean (of means) of CPMs
        if args.library_sizes:
            count_col_idx = 3
            try:
                lib_sizes = [int(open(f).readlines()[0].rstrip()) \
                                    for lib_size_group in args.library_sizes \
                                    for f in lib_size_group]
            except:
                lib_sizes = [int(f) for lib_size_group in args.library_sizes \
                                  for f in lib_size_group]
            counts_tmp.iloc[:, 3:] = counts_tmp.iloc[:, 3:]/lib_sizes*1e6
            means = []
            for lib_size_group in args.library_sizes:
                means.append(counts_tmp.iloc[:, count_col_idx:count_col_idx+len(lib_size_group)]\
                             .mean(axis=1))
                count_col_idx+=len(lib_size_group)
            counts_tmp['accessibility_cpm'] = np.max(means, axis=0)
            print(counts_tmp.head())

        else:
            # assume that the accessibity file has been transformed, rename column
            counts_tmp.columns = counts_tmp.columns.tolist()[:-1] + ["accessibility_cpm"]
        tmp = guide_scan_df.join(counts_tmp.loc[:, [3, "accessibility_cpm"]])
        print(guide_scan_df.columns)
        tmp.drop(3, axis=1, inplace=True)
        guide_scan_df = tmp
        
        
    # Extract peak ID from gRNA label
    guide_scan_df.loc[:, "peak"] = guide_scan_df["gRNA label"].apply(lambda x: "_".join(x.split("_")[:3]))
    
    # Sort guides within peak (by either accessibility or specificity score)
    if args.guides_accesibility_coverage:
        all_grnas_sorted = guide_scan_df.sort_values(by=['peak', 'accessibility_cpm', 'cutting specificity score'],
                                                     ascending=[True, False, False])
    else:
        all_grnas_sorted = guide_scan_df.sort_values(by=['peak','cutting specificity score'],
                                                     ascending=[True, False])

    # Rank gRNAs in peaks controlling the maximum overlap among any two guides
    all_grnas_sorted_nonoverlapping = pd.concat([
        sort_grnas(ii[1], maximum_overlap=args.max_overlap) for ii in all_grnas_sorted.groupby(['peak'])])
    all_grnas_sorted_nonoverlapping['rank_in_peak'] = \
        np.concatenate([range(1, len(ii[1])+1) for ii in all_grnas_sorted_nonoverlapping.groupby(['peak'])])
    all_grnas_sorted_nonoverlapping.to_csv(args.outfile, sep='\t',index=False)


if __name__ == "__main__":
    main()    