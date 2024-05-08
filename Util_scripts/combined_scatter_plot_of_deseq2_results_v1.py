#!/usr/bin/env python
import matplotlib.pyplot as plt
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['lines.markersize'] = 5
import seaborn as sns
import pandas as pd
import numpy as np
from palettable.cartocolors.qualitative import Bold_8_r
import argparse
from itertools import product

def distance_to_genomic_anno(d, promoter_thres=3000):
    if d < promoter_thres:
        return "Promoter <=3kb"
    return "Distal >3kb"

def get_enrich_group(x, fdr_thres=0.01, labels=["DA", "DEG"], padj_cols=["padj_file1", "padj_file2"]):
    if (x[padj_cols[0]] < fdr_thres) and (x[padj_cols[1]] < fdr_thres):
        return 'sig. %(label_file2)s; sig. %(label_file1)s' % ({"label_file1": labels[0], "label_file2": labels[1]})
#         return 'sig. DEG; sig. DA'
    if (x[padj_cols[0]] < fdr_thres) and (x[padj_cols[1]] >= fdr_thres):
        return 'nonsig. %(label_file2)s; sig. %(label_file1)s' % ({"label_file1": labels[0], "label_file2": labels[1]})
#         return 'nonsig. DEG; sig. DA'
    if (x[padj_cols[0]] >= fdr_thres) and (x[padj_cols[1]] < fdr_thres):
        return 'sig. %(label_file2)s; nonsig. %(label_file1)s' % ({"label_file1": labels[0], "label_file2": labels[1]})
#         return 'sig. DEG; nonsig. DA'
    return 'nonsig. %(label_file2)s; nonsig. %(label_file1)s' % ({"label_file1": labels[0], "label_file2": labels[1]})
    
def order_categories_and_colors(categories, colors, fdr_thres=0.01, labels=["DA", "DEG"]):
    cc = []
    cs = []
    for ix, p in enumerate(product(
        ["nonsig. %s" % labels[1], "sig. %s" % labels[1]],
        ["nonsig. %s" % labels[0], "sig. %s" % labels[0]],
    )):
        lab = "; ".join(p)
        if lab in categories:
            cc.append(lab)
            cs.append(colors[ix])
    return cc, cs
    

def main():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter, \
    description="""

    Given two files containing DESeq2 results, create a combined scatter plot
    to compare the results. By default, the script assumes the input files have
    gene ids in the first column, and will join both files accordingly (outer join).

    """)

    ##################################################
    # required args:

    parser.add_argument("-1", "--input-file1", help="""required, file path to first set of DESeq2 results (x-axis).""", 
                        required=True)
    parser.add_argument("-2", "--input-file2", help="""required, file path to second set of DESeq2 results (y-axis).""", 
                        required=True)
    parser.add_argument("-o","--outname", type=str, help="""required, path and basename used to create output files""", 
                        required=True)

    ##################################################
    # optional args:
    parser.add_argument("--colors", nargs='*', default=['grey'] + Bold_8_r.hex_colors[1:2] + [Bold_8_r.hex_colors[3]],
                       help="""optional, list of colors to use.""")
    parser.add_argument("--suffixes", type=str, nargs=2, default=["_file1", "_file2"], 
                        help="""optional, suffixes used to differentiate common columns (default: ['_file1', '_file2'])""")
    parser.add_argument("--padj-cols", type=str, nargs=2, default=["padj_file1", "padj_file2"], 
                        help="""optional, column names of the padj (significance) values (default: ['padj_file1', 'padj_file2'])""")
    parser.add_argument("--log2fc-cols", type=str, nargs=2, default=["log2FoldChange_file1", "log2FoldChange_file2"], 
                        help="""optional, column names of the log2fc (enrichment) values (default: ['log2FoldChange_file1', 'log2FoldChange_file2'])""")
    parser.add_argument("--legend-group-labels", type=str, nargs=2, default=["DA", "DEG"], 
                        help="""optional, tags used to label elements in legend (default: ['DA', 'DEG'])""")
    parser.add_argument("--axes-labels", type=str, nargs=2, default=["DA", "DEG"], 
                        help="""optional, tags used to label axis (default: ['DA', 'DEG'])""")
    parser.add_argument("--suptitle", type=str, default="Combined scatterplot",
                        help="""optional, figure title (default: %(default)s)""")
    parser.add_argument("--how", type=str, default="outer",
                        help="""optional, how to join both result files. Possible options are: {'left', 'right', 'outer', 'inner'},(default: %(default)s)""")
    parser.add_argument("--genes-to-highlight", type=str, nargs="*",
                        help="""optional, gene_symbol names of genes of interest (will be labeled/highlighted)""")
    parser.add_argument("--gene-name-col", type=str, default="GeneName_file1", 
                        help="""optional, column containing the GeneSymbol names used to filter genes to highlight (default: %(default)s)""")
    parser.add_argument("--output-format", type=str, nargs="*", default=["pdf"],
                        help="""optional, extensions used to create plot (e.g. pdf, png). (default: pdf)""")
    parser.add_argument("--fdr-thres", type=float, default=0.05,
                        help="""optional, threshold used to determine significant hits. (default: %(default)f)""")
    parser.add_argument("--input-file1-join-col", type=str, 
                        help="""optional, column name in input file1 used to merge the results. (default: first column)""")
    parser.add_argument("--input-file2-join-col", type=str, 
                        help="""optional, column name in input file2 used to merge the results. (default: first column)""")
    parser.add_argument("--breakdown-by-genomic-anno-distance-col", type=int,
                        help="""optional, if specified, facetwrap the final scatter plot by genomic annotation using the distance
                        specified in this column. This is only useful when the elements from both sets are not the same, 
                        but they are connected (e.g. peaks and genes)""")
    parser.add_argument("--remove-pseudogenes", action='store_true', default=False,
                        help="""optional, the script will filter out "pseudogene"(s) based on the values in the GeneType column (should be present in only one of the files)""")
    parser.add_argument("--fig_size_x", type=float, default=6, help="fig_size_x", dest="fig_size_x")
    parser.add_argument("--fig_size_y", type=float, default=4, help="fig_size_y", dest="fig_size_y")
    
    ###################################################
    
    args = parser.parse_args()
    fdr_thres = args.fdr_thres
    
    # Load files
    df1 = pd.read_csv(args.input_file1, sep='\t')
    df2 = pd.read_csv(args.input_file2, sep='\t')
    
    if args.input_file1_join_col is not None and args.input_file2_join_col is not None:
        merged = df1.merge(df2, 
                           left_on=args.input_file1_join_col, 
                           right_on=args.input_file2_join_col, 
                           suffixes=args.suffixes,
                           how=args.how)
    else:
        merged = df1.join(df2, how=args.how, lsuffix=args.suffixes[0], rsuffix=args.suffixes[1])

    if args.remove_pseudogenes:
        merged = merged.loc[~merged.GeneType.str.contains('pseudogene'), :]
    
    # If there are NA in the significant columns, replace with 1s
#     merged.loc[merged[args.padj_cols[0]].isna(), args.padj_cols[0]] = 1
#     merged.loc[merged[args.padj_cols[1]].isna(), args.padj_cols[1]] = 1
#     merged[args.padj_cols[0]] = merged[args.padj_cols[0]].astype(float)
#     merged[args.padj_cols[1]] = merged[args.padj_cols[1]].astype(float)
    
    # Break down by combination of significant hits
    merged['group (FDR<%.3f)' % fdr_thres] = merged.apply(get_enrich_group, 
                                                          fdr_thres=fdr_thres,
                                                          labels=args.legend_group_labels,
                                                          padj_cols=args.padj_cols,
                                                          axis = 1)

    merged.to_csv('%s.txt' % (args.outname), sep="\t", index=False)
    
    
    corr_coef = merged[args.log2fc_cols[0]].corr(merged[args.log2fc_cols[1]])
               
    hue_order, colors_order = order_categories_and_colors((merged['group (FDR<%.3f)' % fdr_thres]).unique(),                           
                                args.colors, 
                                fdr_thres=args.fdr_thres,
                                labels=args.legend_group_labels)

    fig, ax = plt.subplots(figsize=(args.fig_size_x,args.fig_size_y))
    sns.scatterplot(data = merged, x=args.log2fc_cols[0], y=args.log2fc_cols[1], 
                    hue = 'group (FDR<%.3f)' % fdr_thres, hue_order=hue_order,
                    palette= colors_order, s=20, edgecolor=None, alpha=.75)
    
#     xlims = [np.floor(merged[args.log2fc_cols[0]].min()), np.ceil(merged[args.log2fc_cols[0]].max())]
#     ylims = [np.floor(merged[args.log2fc_cols[1]].min()), np.ceil(merged[args.log2fc_cols[1]].max())]
    xlims = [1.1*(merged[args.log2fc_cols[0]].min()), 1.1*(merged[args.log2fc_cols[0]].max())]
    ylims = [1.1*(merged[args.log2fc_cols[1]].min()), 1.1*(merged[args.log2fc_cols[1]].max())]
    plt.xlim(xlims) 
    plt.ylim(ylims)
    plt.vlines(0, ylims[0], ylims[1], linestyle='--', color='grey')
    plt.hlines(0, xlims[0], xlims[1], linestyle='--', color='grey')

    for x, y, t in zip(
        merged[args.log2fc_cols[0]].values, 
        merged[args.log2fc_cols[1]].values, 
        merged[args.gene_name_col]) :
        if t not in args.genes_to_highlight: continue
        ax.annotate('{}'.format(t), xy=(x, y), 
                    xytext=(.75*x, .75*y), #xytext=(5, 0), 
                    ha='left', va='bottom', #ha='left',
                    arrowprops=dict(arrowstyle='-', color='black'))

    plt.legend(bbox_to_anchor=(1.05, 1), frameon=False, fontsize=8)
    plt.ylabel(args.axes_labels[1])#'RNA-seq log2fc')
    plt.xlabel(args.axes_labels[0])#'K27ac ChIP-seq log2fc')
    fig.suptitle('%s for log2FC values; correlation coef = %.3f' % (args.suptitle, corr_coef), fontsize=12)#'Cre+/g5 v Cre+/NTC')

    sns.despine()
    plt.tight_layout()
    for extension in args.output_format:
        plt.savefig('%s.%s' % (args.outname, extension))
    
    return

if __name__ == '__main__':
    main()
