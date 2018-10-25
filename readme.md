# Cooperative TMH analysis source code.

## Running the code

### Requirements

The script requires Python3, SciPy, biopython, Numpy, matplotlib, git, and git-lfs.
TMSOC scripts provided by Dr Wing Cheong Wong. The software is explained and demonstrated at [tmsoc.bii.a-star.edu.sg](http://tmsoc.bii.a-star.edu.sg/).

(Wong, W.-C., Maurer-Stroh, S., Schneider, G., & Eisenhaber, F. (2012). Transmembrane helix: simple or complex. Nucleic Acids Research, 40(W1), W370–W375. <https://doi.org/10.1093/nar/gks379>)

### Generating the graphs

1. Download the files using `git clone https://github.com/JamesABaker/multipass_tmh_analysis.git`.
2. From within `multipass_tmh_analysis`, run `git lfs install` and `git lfs fetch`.
3. Move your uniprot text file into the `scripts/` folder and then run:

- Usage

      python3 [Script] [file] [TMH number for graph axis] [OPTIONAL TMH filter number]

- Example

      python3 Uniprot_to_graphs_and_stats.py file.txt 7 7

Please log any issues you run into on the [GitHub issue page](https://github.com/JamesABaker/multipass_tmh_analysis/issues/new).

## File structure

The UniProt txt file should be moved to the current working directory as well as all the files from `scripts/`.

Raw images from the original output are also included, as well as the original zip file downloads from UniProt which often tweaks the database.

    .
    ├── KW-zip-files-uniref50-swissprot
    │   ├── Integrin.txt.gz
    │   ├── Ion\ channel.txt.gz
    │   ├── Potassium\ channel.gz
    │   ├── calcium\ channel.txt.gz
    │   ├── chloride\ channel.gz
    │   ├── ion\ transport.txt.gz
    │   ├── ligand\ gated\ ion\ channel.txt.gz
    │   ├── porin.txt.gz
    │   ├── protein\ transport.txt.gz
    │   ├── sodium\ channel.txt.gz
    │   ├── sugar\ transport.txt.gz
    │   ├── transmembrane\ keywords\ query\ results.md
    │   ├── uniref-uniprot%3A%28keyword%3A%22Calcium+channel+%5BKW-0107%5D%22+annotation%3A%--.list
    │   ├── uniref-uniprot%3A%28keyword%3A%22Calcium+channel+%5BKW-0107%5D%22+annotation%3A%--.list.gz
    │   ├── uniref-uniprot%3A%28keyword%3A%22Chloride+channel+%5BKW-0869%5D%22+annotation%3A--.list
    │   ├── uniref-uniprot%3A%28keyword%3A%22Chloride+channel+%5BKW-0869%5D%22+annotation%3A--.list.gz
    │   ├── uniref-uniprot%3A%28keyword%3A%22Integrin+%5BKW-0401%5D%22+annotation%3A%28type%--.list
    │   ├── uniref-uniprot%3A%28keyword%3A%22Integrin+%5BKW-0401%5D%22+annotation%3A%28type%--.list.gz
    │   ├── uniref-uniprot%3A%28keyword%3A%22Ion+channel+%5BKW-0407%5D%22+annotation%3A%28ty--.list
    │   ├── uniref-uniprot%3A%28keyword%3A%22Ion+channel+%5BKW-0407%5D%22+annotation%3A%28ty--.list.gz
    │   ├── uniref-uniprot%3A%28keyword%3A%22Ion+transport+%5BKW-0406%5D%22+annotation%3A%28--.list
    │   ├── uniref-uniprot%3A%28keyword%3A%22Ion+transport+%5BKW-0406%5D%22+annotation%3A%28--.list.gz
    │   ├── uniref-uniprot%3A%28keyword%3A%22Ligand-gated+ion+channel+%5BKW-1071%5D%22+annot--.list
    │   ├── uniref-uniprot%3A%28keyword%3A%22Ligand-gated+ion+channel+%5BKW-1071%5D%22+annot--.list.gz
    │   ├── uniref-uniprot%3A%28keyword%3A%22Notch+signaling+pathway+%5BKW-0914%5D%22+annota--.list
    │   ├── uniref-uniprot%3A%28keyword%3A%22Notch+signaling+pathway+%5BKW-0914%5D%22+annota--.list.gz
    │   ├── uniref-uniprot%3A%28keyword%3A%22Porin+%5BKW-0626%5D%22+annotation%3A%28type%3At--.list
    │   ├── uniref-uniprot%3A%28keyword%3A%22Porin+%5BKW-0626%5D%22+annotation%3A%28type%3At--.list\ 16.34.16.gz
    │   ├── uniref-uniprot%3A%28keyword%3A%22Potassium+channel+%5BKW-0631%5D%22+annotation%3--.list
    │   ├── uniref-uniprot%3A%28keyword%3A%22Potassium+channel+%5BKW-0631%5D%22+annotation%3--.list.gz
    │   ├── uniref-uniprot%3A%28keyword%3A%22Protein+transport+%5BKW-0653%5D%22+annotation%3--\ 2.list
    │   ├── uniref-uniprot%3A%28keyword%3A%22Protein+transport+%5BKW-0653%5D%22+annotation%3--.list.gz
    │   ├── uniref-uniprot%3A%28keyword%3A%22Sodium+channel+%5BKW-0894%5D%22+annotation%3A%2--.list
    │   ├── uniref-uniprot%3A%28keyword%3A%22Sodium+channel+%5BKW-0894%5D%22+annotation%3A%2--.list.gz
    │   ├── uniref-uniprot%3A%28keyword%3A%22Sugar+transport+%5BKW-0762%5D%22+annotation%3A%--\ 2.list
    │   ├── uniref-uniprot%3A%28keyword%3A%22Sugar+transport+%5BKW-0762%5D%22+annotation%3A%--.list
    │   ├── uniref-uniprot%3A%28keyword%3A%22Sugar+transport+%5BKW-0762%5D%22+annotation%3A%--.list.gz
    │   ├── uniref-uniprot%3A%28keyword%3A%22Voltage-gated+channel+%5BKW-0851%5D%22+annotati--.list
    │   ├── uniref-uniprot%3A%28keyword%3A%22Voltage-gated+channel+%5BKW-0851%5D%22+annotati--.list.gz
    │   └── voltage\ gated\ ion\ channel.txt.gz
    ├── images
    │   ├── 00_47_46__21_12_2017_GPCR_UniRef50_test.txt_Complexity.pdf
    │   ├── 00_47_46__21_12_2017_GPCR_UniRef50_test.txt_Complexity.png
    │   ├── 00_47_47__21_12_2017_GPCR_UniRef50_test.txt_Length.pdf
    │   ├── 00_47_47__21_12_2017_GPCR_UniRef50_test.txt_Length.png
    │   ├── 00_47_48__21_12_2017_GPCR_UniRef50_test.txt_Eisenberg\ hydrophobiciity\ scale.pdf
    │   ├── 00_47_48__21_12_2017_GPCR_UniRef50_test.txt_Eisenberg\ hydrophobiciity\ scale.png
    │   ├── 00_47_48__21_12_2017_GPCR_UniRef50_test.txt_Kyte\ &\ Doolittle\ hydrophobiciity\ scale.pdf
    │   ├── 00_47_48__21_12_2017_GPCR_UniRef50_test.txt_Kyte\ &\ Doolittle\ hydrophobiciity\ scale.png
    │   ├── 00_47_49__21_12_2017_GPCR_UniRef50_test.txt_\ White\ and\ Wimley\ hydrophobiciity\ scale.pdf
    │   ├── 00_47_49__21_12_2017_GPCR_UniRef50_test.txt_\ White\ and\ Wimley\ hydrophobiciity\ scale.png
    │   ├── 09_38_19__02_02_2018_KW\ potassium\ channel.txt_Complexity.pdf
    │   ├── 09_38_19__02_02_2018_KW\ potassium\ channel.txt_Complexity.png
    │   ├── 09_38_21__02_02_2018_KW\ potassium\ channel.txt_Length.pdf
    │   ├── 09_38_21__02_02_2018_KW\ potassium\ channel.txt_Length.png
    │   ├── 09_38_23__02_02_2018_KW\ potassium\ channel.txt_Kyte\ &\ Doolittle\ hydrophobiciity\ scale.pdf
    │   ├── 09_38_23__02_02_2018_KW\ potassium\ channel.txt_Kyte\ &\ Doolittle\ hydrophobiciity\ scale.png
    │   ├── 09_38_25__02_02_2018_KW\ potassium\ channel.txt_Eisenberg\ hydrophobiciity\ scale.pdf
    │   ├── 09_38_25__02_02_2018_KW\ potassium\ channel.txt_Eisenberg\ hydrophobiciity\ scale.png
    │   ├── 09_38_27__02_02_2018_KW\ potassium\ channel.txt_\ White\ and\ Wimley\ hydrophobiciity\ scale.pdf
    │   ├── 09_38_27__02_02_2018_KW\ potassium\ channel.txt_\ White\ and\ Wimley\ hydrophobiciity\ scale.png
    │   ├── 09_48_37__02_02_2018_KW\ calcium\ channel.txt_Complexity.pdf
    │   ├── 09_48_37__02_02_2018_KW\ calcium\ channel.txt_Complexity.png
    │   ├── 09_48_39__02_02_2018_KW\ calcium\ channel.txt_Length.pdf
    │   ├── 09_48_39__02_02_2018_KW\ calcium\ channel.txt_Length.png
    │   ├── 09_48_41__02_02_2018_KW\ calcium\ channel.txt_Kyte\ &\ Doolittle\ hydrophobiciity\ scale.pdf
    │   ├── 09_48_41__02_02_2018_KW\ calcium\ channel.txt_Kyte\ &\ Doolittle\ hydrophobiciity\ scale.png
    │   ├── 09_48_43__02_02_2018_KW\ calcium\ channel.txt_Eisenberg\ hydrophobiciity\ scale.pdf
    │   ├── 09_48_43__02_02_2018_KW\ calcium\ channel.txt_Eisenberg\ hydrophobiciity\ scale.png
    │   ├── 09_48_45__02_02_2018_KW\ calcium\ channel.txt_\ White\ and\ Wimley\ hydrophobiciity\ scale.pdf
    │   ├── 09_48_45__02_02_2018_KW\ calcium\ channel.txt_\ White\ and\ Wimley\ hydrophobiciity\ scale.png
    │   ├── 10_00_05__02_02_2018_KW\ chloride\ channel.txt_Complexity.pdf
    │   ├── 10_00_05__02_02_2018_KW\ chloride\ channel.txt_Complexity.png
    │   ├── 10_00_07__02_02_2018_KW\ chloride\ channel.txt_Length.pdf
    │   ├── 10_00_07__02_02_2018_KW\ chloride\ channel.txt_Length.png
    │   ├── 10_00_09__02_02_2018_KW\ chloride\ channel.txt_Kyte\ &\ Doolittle\ hydrophobiciity\ scale.pdf
    │   ├── 10_00_09__02_02_2018_KW\ chloride\ channel.txt_Kyte\ &\ Doolittle\ hydrophobiciity\ scale.png
    │   ├── 10_00_11__02_02_2018_KW\ chloride\ channel.txt_Eisenberg\ hydrophobiciity\ scale.pdf
    │   ├── 10_00_11__02_02_2018_KW\ chloride\ channel.txt_Eisenberg\ hydrophobiciity\ scale.png
    │   ├── 10_00_13__02_02_2018_KW\ chloride\ channel.txt_\ White\ and\ Wimley\ hydrophobiciity\ scale.pdf
    │   ├── 10_00_13__02_02_2018_KW\ chloride\ channel.txt_\ White\ and\ Wimley\ hydrophobiciity\ scale.png
    │   ├── 10_03_33__02_02_2018_KW\ calcium\ channel.txt_Complexity.pdf
    │   ├── 10_03_33__02_02_2018_KW\ calcium\ channel.txt_Complexity.png
    │   ├── 10_03_35__02_02_2018_KW\ calcium\ channel.txt_Length.pdf
    │   ├── 10_03_35__02_02_2018_KW\ calcium\ channel.txt_Length.png
    │   ├── 10_03_37__02_02_2018_KW\ calcium\ channel.txt_Kyte\ &\ Doolittle\ hydrophobiciity\ scale.pdf
    │   ├── 10_03_37__02_02_2018_KW\ calcium\ channel.txt_Kyte\ &\ Doolittle\ hydrophobiciity\ scale.png
    │   ├── 10_03_39__02_02_2018_KW\ calcium\ channel.txt_Eisenberg\ hydrophobiciity\ scale.pdf
    │   ├── 10_03_39__02_02_2018_KW\ calcium\ channel.txt_Eisenberg\ hydrophobiciity\ scale.png
    │   ├── 10_03_41__02_02_2018_KW\ calcium\ channel.txt_\ White\ and\ Wimley\ hydrophobiciity\ scale.pdf
    │   ├── 10_03_41__02_02_2018_KW\ calcium\ channel.txt_\ White\ and\ Wimley\ hydrophobiciity\ scale.png
    │   ├── 10_14_43__02_02_2018_KW\ calcium\ channel.txt_Complexity.pdf
    │   ├── 10_14_43__02_02_2018_KW\ calcium\ channel.txt_Complexity.png
    │   ├── 10_14_45__02_02_2018_KW\ calcium\ channel.txt_Length.pdf
    │   ├── 10_14_45__02_02_2018_KW\ calcium\ channel.txt_Length.png
    │   ├── 10_14_47__02_02_2018_KW\ calcium\ channel.txt_Kyte\ &\ Doolittle\ hydrophobiciity\ scale.pdf
    │   ├── 10_14_47__02_02_2018_KW\ calcium\ channel.txt_Kyte\ &\ Doolittle\ hydrophobiciity\ scale.png
    │   ├── 10_14_49__02_02_2018_KW\ calcium\ channel.txt_Eisenberg\ hydrophobiciity\ scale.pdf
    │   ├── 10_14_49__02_02_2018_KW\ calcium\ channel.txt_Eisenberg\ hydrophobiciity\ scale.png
    │   ├── 10_14_51__02_02_2018_KW\ calcium\ channel.txt_\ White\ and\ Wimley\ hydrophobiciity\ scale.pdf
    │   ├── 10_14_51__02_02_2018_KW\ calcium\ channel.txt_\ White\ and\ Wimley\ hydrophobiciity\ scale.png
    │   ├── 10_24_28__02_02_2018_KW\ chloride\ channel.txt_Complexity.pdf
    │   ├── 10_24_28__02_02_2018_KW\ chloride\ channel.txt_Complexity.png
    │   ├── 10_24_30__02_02_2018_KW\ chloride\ channel.txt_Length.pdf
    │   ├── 10_24_30__02_02_2018_KW\ chloride\ channel.txt_Length.png
    │   ├── 10_24_32__02_02_2018_KW\ chloride\ channel.txt_Kyte\ &\ Doolittle\ hydrophobiciity\ scale.pdf
    │   ├── 10_24_32__02_02_2018_KW\ chloride\ channel.txt_Kyte\ &\ Doolittle\ hydrophobiciity\ scale.png
    │   ├── 10_24_34__02_02_2018_KW\ chloride\ channel.txt_Eisenberg\ hydrophobiciity\ scale.pdf
    │   ├── 10_24_34__02_02_2018_KW\ chloride\ channel.txt_Eisenberg\ hydrophobiciity\ scale.png
    │   ├── 10_24_36__02_02_2018_KW\ chloride\ channel.txt_\ White\ and\ Wimley\ hydrophobiciity\ scale.pdf
    │   ├── 10_24_36__02_02_2018_KW\ chloride\ channel.txt_\ White\ and\ Wimley\ hydrophobiciity\ scale.png
    │   ├── 13_00_37__02_02_2018_KW\ sugar\ transport.txt_Complexity.pdf
    │   ├── 13_00_37__02_02_2018_KW\ sugar\ transport.txt_Complexity.png
    │   ├── 13_00_39__02_02_2018_KW\ sugar\ transport.txt_Length.pdf
    │   ├── 13_00_39__02_02_2018_KW\ sugar\ transport.txt_Length.png
    │   ├── 13_00_41__02_02_2018_KW\ sugar\ transport.txt_Kyte\ &\ Doolittle\ hydrophobiciity\ scale.pdf
    │   ├── 13_00_41__02_02_2018_KW\ sugar\ transport.txt_Kyte\ &\ Doolittle\ hydrophobiciity\ scale.png
    │   ├── 13_00_43__02_02_2018_KW\ sugar\ transport.txt_Eisenberg\ hydrophobiciity\ scale.pdf
    │   ├── 13_00_43__02_02_2018_KW\ sugar\ transport.txt_Eisenberg\ hydrophobiciity\ scale.png
    │   ├── 13_00_45__02_02_2018_KW\ sugar\ transport.txt_\ White\ and\ Wimley\ hydrophobiciity\ scale.pdf
    │   ├── 13_00_45__02_02_2018_KW\ sugar\ transport.txt_\ White\ and\ Wimley\ hydrophobiciity\ scale.png
    │   ├── 13_01_46__02_02_2018_KW\ sugar\ transport.txt_Complexity.pdf
    │   ├── 13_01_46__02_02_2018_KW\ sugar\ transport.txt_Complexity.png
    │   ├── 13_01_48__02_02_2018_KW\ sugar\ transport.txt_Length.pdf
    │   ├── 13_01_48__02_02_2018_KW\ sugar\ transport.txt_Length.png
    │   ├── 13_01_49__02_02_2018_KW\ sugar\ transport.txt_Kyte\ &\ Doolittle\ hydrophobiciity\ scale.pdf
    │   ├── 13_01_49__02_02_2018_KW\ sugar\ transport.txt_Kyte\ &\ Doolittle\ hydrophobiciity\ scale.png
    │   ├── 13_01_51__02_02_2018_KW\ sugar\ transport.txt_Eisenberg\ hydrophobiciity\ scale.pdf
    │   ├── 13_01_51__02_02_2018_KW\ sugar\ transport.txt_Eisenberg\ hydrophobiciity\ scale.png
    │   ├── 13_01_53__02_02_2018_KW\ sugar\ transport.txt_\ White\ and\ Wimley\ hydrophobiciity\ scale.pdf
    │   ├── 13_01_53__02_02_2018_KW\ sugar\ transport.txt_\ White\ and\ Wimley\ hydrophobiciity\ scale.png
    │   ├── 13_03_03__02_02_2018_KW\ sugar\ transport.txt_Complexity.pdf
    │   ├── 13_03_03__02_02_2018_KW\ sugar\ transport.txt_Complexity.png
    │   ├── 13_03_05__02_02_2018_KW\ sugar\ transport.txt_Length.pdf
    │   ├── 13_03_05__02_02_2018_KW\ sugar\ transport.txt_Length.png
    │   ├── 13_03_07__02_02_2018_KW\ sugar\ transport.txt_Kyte\ &\ Doolittle\ hydrophobiciity\ scale.pdf
    │   ├── 13_03_07__02_02_2018_KW\ sugar\ transport.txt_Kyte\ &\ Doolittle\ hydrophobiciity\ scale.png
    │   ├── 13_03_09__02_02_2018_KW\ sugar\ transport.txt_Eisenberg\ hydrophobiciity\ scale.pdf
    │   ├── 13_03_09__02_02_2018_KW\ sugar\ transport.txt_Eisenberg\ hydrophobiciity\ scale.png
    │   ├── 13_03_11__02_02_2018_KW\ sugar\ transport.txt_\ White\ and\ Wimley\ hydrophobiciity\ scale.pdf
    │   ├── 13_03_11__02_02_2018_KW\ sugar\ transport.txt_\ White\ and\ Wimley\ hydrophobiciity\ scale.png
    │   ├── 13_08_45__31_01_2018_KW\ porin.txt_Complexity.pdf
    │   ├── 13_08_45__31_01_2018_KW\ porin.txt_Complexity.png
    │   ├── 13_08_47__31_01_2018_KW\ porin.txt_Length.pdf
    │   ├── 13_08_47__31_01_2018_KW\ porin.txt_Length.png
    │   ├── 13_08_49__31_01_2018_KW\ porin.txt_Kyte\ &\ Doolittle\ hydrophobiciity\ scale.pdf
    │   ├── 13_08_49__31_01_2018_KW\ porin.txt_Kyte\ &\ Doolittle\ hydrophobiciity\ scale.png
    │   ├── 13_08_50__31_01_2018_KW\ porin.txt_Eisenberg\ hydrophobiciity\ scale.pdf
    │   ├── 13_08_50__31_01_2018_KW\ porin.txt_Eisenberg\ hydrophobiciity\ scale.png
    │   ├── 13_08_52__31_01_2018_KW\ porin.txt_\ White\ and\ Wimley\ hydrophobiciity\ scale.pdf
    │   ├── 13_08_52__31_01_2018_KW\ porin.txt_\ White\ and\ Wimley\ hydrophobiciity\ scale.png
    │   ├── 13_27_58__31_01_2018_KW\ potassium\ channel.txt_Complexity.pdf
    │   ├── 13_27_58__31_01_2018_KW\ potassium\ channel.txt_Complexity.png
    │   ├── 13_28_00__31_01_2018_KW\ potassium\ channel.txt_Length.pdf
    │   ├── 13_28_00__31_01_2018_KW\ potassium\ channel.txt_Length.png
    │   ├── 13_28_02__31_01_2018_KW\ potassium\ channel.txt_Kyte\ &\ Doolittle\ hydrophobiciity\ scale.pdf
    │   ├── 13_28_02__31_01_2018_KW\ potassium\ channel.txt_Kyte\ &\ Doolittle\ hydrophobiciity\ scale.png
    │   ├── 13_28_03__31_01_2018_KW\ potassium\ channel.txt_Eisenberg\ hydrophobiciity\ scale.pdf
    │   ├── 13_28_03__31_01_2018_KW\ potassium\ channel.txt_Eisenberg\ hydrophobiciity\ scale.png
    │   ├── 13_28_05__31_01_2018_KW\ potassium\ channel.txt_\ White\ and\ Wimley\ hydrophobiciity\ scale.pdf
    │   ├── 13_28_05__31_01_2018_KW\ potassium\ channel.txt_\ White\ and\ Wimley\ hydrophobiciity\ scale.png
    │   ├── 13_32_28__31_01_2018_KW\ ion\ channel.txt_Complexity.pdf
    │   ├── 13_32_28__31_01_2018_KW\ ion\ channel.txt_Complexity.png
    │   ├── 13_32_30__31_01_2018_KW\ ion\ channel.txt_Length.pdf
    │   ├── 13_32_30__31_01_2018_KW\ ion\ channel.txt_Length.png
    │   ├── 13_32_32__31_01_2018_KW\ ion\ channel.txt_Kyte\ &\ Doolittle\ hydrophobiciity\ scale.pdf
    │   ├── 13_32_32__31_01_2018_KW\ ion\ channel.txt_Kyte\ &\ Doolittle\ hydrophobiciity\ scale.png
    │   ├── 13_32_34__31_01_2018_KW\ ion\ channel.txt_Eisenberg\ hydrophobiciity\ scale.pdf
    │   ├── 13_32_34__31_01_2018_KW\ ion\ channel.txt_Eisenberg\ hydrophobiciity\ scale.png
    │   ├── 13_32_35__31_01_2018_KW\ ion\ channel.txt_\ White\ and\ Wimley\ hydrophobiciity\ scale.pdf
    │   ├── 13_32_35__31_01_2018_KW\ ion\ channel.txt_\ White\ and\ Wimley\ hydrophobiciity\ scale.png
    │   ├── 13_40_22__31_01_2018_KW\ potassium\ channel.txt_Complexity.pdf
    │   ├── 13_40_22__31_01_2018_KW\ potassium\ channel.txt_Complexity.png
    │   ├── 13_40_24__31_01_2018_KW\ potassium\ channel.txt_Length.pdf
    │   ├── 13_40_24__31_01_2018_KW\ potassium\ channel.txt_Length.png
    │   ├── 13_40_26__31_01_2018_KW\ potassium\ channel.txt_Kyte\ &\ Doolittle\ hydrophobiciity\ scale.pdf
    │   ├── 13_40_26__31_01_2018_KW\ potassium\ channel.txt_Kyte\ &\ Doolittle\ hydrophobiciity\ scale.png
    │   ├── 13_40_28__31_01_2018_KW\ potassium\ channel.txt_Eisenberg\ hydrophobiciity\ scale.pdf
    │   ├── 13_40_28__31_01_2018_KW\ potassium\ channel.txt_Eisenberg\ hydrophobiciity\ scale.png
    │   ├── 13_40_29__31_01_2018_KW\ potassium\ channel.txt_\ White\ and\ Wimley\ hydrophobiciity\ scale.pdf
    │   ├── 13_40_29__31_01_2018_KW\ potassium\ channel.txt_\ White\ and\ Wimley\ hydrophobiciity\ scale.png
    │   ├── 14_04_02__25_01_2018_KW\ ion\ channel.txt_Complexity.pdf
    │   ├── 14_04_02__25_01_2018_KW\ ion\ channel.txt_Complexity.png
    │   ├── 14_04_03__25_01_2018_KW\ ion\ channel.txt_Length.pdf
    │   ├── 14_04_03__25_01_2018_KW\ ion\ channel.txt_Length.png
    │   ├── 14_04_04__25_01_2018_KW\ ion\ channel.txt_Kyte\ &\ Doolittle\ hydrophobiciity\ scale.pdf
    │   ├── 14_04_04__25_01_2018_KW\ ion\ channel.txt_Kyte\ &\ Doolittle\ hydrophobiciity\ scale.png
    │   ├── 14_04_05__25_01_2018_KW\ ion\ channel.txt_\ White\ and\ Wimley\ hydrophobiciity\ scale.pdf
    │   ├── 14_04_05__25_01_2018_KW\ ion\ channel.txt_\ White\ and\ Wimley\ hydrophobiciity\ scale.png
    │   ├── 14_04_05__25_01_2018_KW\ ion\ channel.txt_Eisenberg\ hydrophobiciity\ scale.pdf
    │   ├── 14_04_05__25_01_2018_KW\ ion\ channel.txt_Eisenberg\ hydrophobiciity\ scale.png
    │   ├── 14_07_57__25_01_2018_7TM_list_UniRef50.txt_Complexity.pdf
    │   ├── 14_07_57__25_01_2018_7TM_list_UniRef50.txt_Complexity.png
    │   ├── 14_07_58__25_01_2018_7TM_list_UniRef50.txt_Length.pdf
    │   ├── 14_07_58__25_01_2018_7TM_list_UniRef50.txt_Length.png
    │   ├── 14_07_59__25_01_2018_7TM_list_UniRef50.txt_Eisenberg\ hydrophobiciity\ scale.pdf
    │   ├── 14_07_59__25_01_2018_7TM_list_UniRef50.txt_Eisenberg\ hydrophobiciity\ scale.png
    │   ├── 14_07_59__25_01_2018_7TM_list_UniRef50.txt_Kyte\ &\ Doolittle\ hydrophobiciity\ scale.pdf
    │   ├── 14_07_59__25_01_2018_7TM_list_UniRef50.txt_Kyte\ &\ Doolittle\ hydrophobiciity\ scale.png
    │   ├── 14_08_00__25_01_2018_7TM_list_UniRef50.txt_\ White\ and\ Wimley\ hydrophobiciity\ scale.pdf
    │   ├── 14_08_00__25_01_2018_7TM_list_UniRef50.txt_\ White\ and\ Wimley\ hydrophobiciity\ scale.png
    │   ├── 17_16_40__24_01_2018_voltage_gated_ion_channel_uniref50.txt_Complexity.pdf
    │   ├── 17_16_40__24_01_2018_voltage_gated_ion_channel_uniref50.txt_Complexity.png
    │   ├── 17_16_42__24_01_2018_voltage_gated_ion_channel_uniref50.txt_Length.pdf
    │   ├── 17_16_42__24_01_2018_voltage_gated_ion_channel_uniref50.txt_Length.png
    │   ├── 17_16_44__24_01_2018_voltage_gated_ion_channel_uniref50.txt_Kyte\ &\ Doolittle\ hydrophobiciity\ scale.pdf
    │   ├── 17_16_44__24_01_2018_voltage_gated_ion_channel_uniref50.txt_Kyte\ &\ Doolittle\ hydrophobiciity\ scale.png
    │   ├── 17_16_46__24_01_2018_voltage_gated_ion_channel_uniref50.txt_Eisenberg\ hydrophobiciity\ scale.pdf
    │   ├── 17_16_46__24_01_2018_voltage_gated_ion_channel_uniref50.txt_Eisenberg\ hydrophobiciity\ scale.png
    │   ├── 17_16_48__24_01_2018_voltage_gated_ion_channel_uniref50.txt_\ White\ and\ Wimley\ hydrophobiciity\ scale.pdf
    │   ├── 17_16_48__24_01_2018_voltage_gated_ion_channel_uniref50.txt_\ White\ and\ Wimley\ hydrophobiciity\ scale.png
    │   ├── 17_33_11__30_01_2018_KW\ porin.txt_Complexity.pdf
    │   ├── 17_33_11__30_01_2018_KW\ porin.txt_Complexity.png
    │   ├── 17_33_12__30_01_2018_KW\ porin.txt_Length.pdf
    │   ├── 17_33_12__30_01_2018_KW\ porin.txt_Length.png
    │   ├── 17_33_13__30_01_2018_KW\ porin.txt_Eisenberg\ hydrophobiciity\ scale.pdf
    │   ├── 17_33_13__30_01_2018_KW\ porin.txt_Eisenberg\ hydrophobiciity\ scale.png
    │   ├── 17_33_13__30_01_2018_KW\ porin.txt_Kyte\ &\ Doolittle\ hydrophobiciity\ scale.pdf
    │   ├── 17_33_13__30_01_2018_KW\ porin.txt_Kyte\ &\ Doolittle\ hydrophobiciity\ scale.png
    │   ├── 17_33_14__30_01_2018_KW\ porin.txt_\ White\ and\ Wimley\ hydrophobiciity\ scale.pdf
    │   ├── 17_33_14__30_01_2018_KW\ porin.txt_\ White\ and\ Wimley\ hydrophobiciity\ scale.png
    │   ├── 17_36_28__30_01_2018_KW\ porin.txt_Complexity.pdf
    │   ├── 17_36_28__30_01_2018_KW\ porin.txt_Complexity.png
    │   ├── 17_36_29__30_01_2018_KW\ porin.txt_Length.pdf
    │   ├── 17_36_29__30_01_2018_KW\ porin.txt_Length.png
    │   ├── 17_36_30__30_01_2018_KW\ porin.txt_Kyte\ &\ Doolittle\ hydrophobiciity\ scale.pdf
    │   ├── 17_36_30__30_01_2018_KW\ porin.txt_Kyte\ &\ Doolittle\ hydrophobiciity\ scale.png
    │   ├── 17_36_31__30_01_2018_KW\ porin.txt_\ White\ and\ Wimley\ hydrophobiciity\ scale.pdf
    │   ├── 17_36_31__30_01_2018_KW\ porin.txt_\ White\ and\ Wimley\ hydrophobiciity\ scale.png
    │   ├── 17_36_31__30_01_2018_KW\ porin.txt_Eisenberg\ hydrophobiciity\ scale.pdf
    │   ├── 17_36_31__30_01_2018_KW\ porin.txt_Eisenberg\ hydrophobiciity\ scale.png
    │   ├── 17_36_42__30_01_2018_KW\ porin.txt_Complexity.pdf
    │   ├── 17_36_42__30_01_2018_KW\ porin.txt_Complexity.png
    │   ├── 17_36_43__30_01_2018_KW\ porin.txt_Length.pdf
    │   ├── 17_36_43__30_01_2018_KW\ porin.txt_Length.png
    │   ├── 17_36_44__30_01_2018_KW\ porin.txt_Kyte\ &\ Doolittle\ hydrophobiciity\ scale.pdf
    │   ├── 17_36_44__30_01_2018_KW\ porin.txt_Kyte\ &\ Doolittle\ hydrophobiciity\ scale.png
    │   ├── 17_36_45__30_01_2018_KW\ porin.txt_Eisenberg\ hydrophobiciity\ scale.pdf
    │   ├── 17_36_45__30_01_2018_KW\ porin.txt_Eisenberg\ hydrophobiciity\ scale.png
    │   ├── 17_36_46__30_01_2018_KW\ porin.txt_\ White\ and\ Wimley\ hydrophobiciity\ scale.pdf
    │   ├── 17_36_46__30_01_2018_KW\ porin.txt_\ White\ and\ Wimley\ hydrophobiciity\ scale.png
    │   ├── 17_37_52__30_01_2018_KW\ porin.txt_Complexity.pdf
    │   ├── 17_37_52__30_01_2018_KW\ porin.txt_Complexity.png
    │   ├── 17_37_52__30_01_2018_KW\ porin.txt_Length.pdf
    │   ├── 17_37_52__30_01_2018_KW\ porin.txt_Length.png
    │   ├── 17_37_53__30_01_2018_KW\ porin.txt_Kyte\ &\ Doolittle\ hydrophobiciity\ scale.pdf
    │   ├── 17_37_53__30_01_2018_KW\ porin.txt_Kyte\ &\ Doolittle\ hydrophobiciity\ scale.png
    │   ├── 17_37_54__30_01_2018_KW\ porin.txt_Eisenberg\ hydrophobiciity\ scale.pdf
    │   ├── 17_37_54__30_01_2018_KW\ porin.txt_Eisenberg\ hydrophobiciity\ scale.png
    │   ├── 17_37_55__30_01_2018_KW\ porin.txt_\ White\ and\ Wimley\ hydrophobiciity\ scale.pdf
    │   ├── 17_37_55__30_01_2018_KW\ porin.txt_\ White\ and\ Wimley\ hydrophobiciity\ scale.png
    │   ├── 17_38_41__30_01_2018_KW\ porin.txt_Complexity.pdf
    │   ├── 17_38_41__30_01_2018_KW\ porin.txt_Complexity.png
    │   ├── 17_38_42__30_01_2018_KW\ porin.txt_Length.pdf
    │   ├── 17_38_42__30_01_2018_KW\ porin.txt_Length.png
    │   ├── 17_38_43__30_01_2018_KW\ porin.txt_Kyte\ &\ Doolittle\ hydrophobiciity\ scale.pdf
    │   ├── 17_38_43__30_01_2018_KW\ porin.txt_Kyte\ &\ Doolittle\ hydrophobiciity\ scale.png
    │   ├── 17_38_44__30_01_2018_KW\ porin.txt_\ White\ and\ Wimley\ hydrophobiciity\ scale.pdf
    │   ├── 17_38_44__30_01_2018_KW\ porin.txt_\ White\ and\ Wimley\ hydrophobiciity\ scale.png
    │   ├── 17_38_44__30_01_2018_KW\ porin.txt_Eisenberg\ hydrophobiciity\ scale.pdf
    │   ├── 17_38_44__30_01_2018_KW\ porin.txt_Eisenberg\ hydrophobiciity\ scale.png
    │   ├── 17_38_55__30_01_2018_KW\ porin.txt_Complexity.pdf
    │   ├── 17_38_55__30_01_2018_KW\ porin.txt_Complexity.png
    │   ├── 17_38_56__30_01_2018_KW\ porin.txt_Length.pdf
    │   ├── 17_38_56__30_01_2018_KW\ porin.txt_Length.png
    │   ├── 17_38_57__30_01_2018_KW\ porin.txt_Kyte\ &\ Doolittle\ hydrophobiciity\ scale.pdf
    │   ├── 17_38_57__30_01_2018_KW\ porin.txt_Kyte\ &\ Doolittle\ hydrophobiciity\ scale.png
    │   ├── 17_38_58__30_01_2018_KW\ porin.txt_Eisenberg\ hydrophobiciity\ scale.pdf
    │   ├── 17_38_58__30_01_2018_KW\ porin.txt_Eisenberg\ hydrophobiciity\ scale.png
    │   ├── 17_38_59__30_01_2018_KW\ porin.txt_\ White\ and\ Wimley\ hydrophobiciity\ scale.pdf
    │   ├── 17_38_59__30_01_2018_KW\ porin.txt_\ White\ and\ Wimley\ hydrophobiciity\ scale.png
    │   ├── 17_43_29__24_01_2018_voltage_gated_ion_channel_uniref50.txt_Complexity.pdf
    │   ├── 17_43_29__24_01_2018_voltage_gated_ion_channel_uniref50.txt_Complexity.png
    │   ├── 17_43_32__24_01_2018_voltage_gated_ion_channel_uniref50.txt_Length.pdf
    │   ├── 17_43_32__24_01_2018_voltage_gated_ion_channel_uniref50.txt_Length.png
    │   ├── 17_43_35__24_01_2018_voltage_gated_ion_channel_uniref50.txt_Kyte\ &\ Doolittle\ hydrophobiciity\ scale.pdf
    │   ├── 17_43_35__24_01_2018_voltage_gated_ion_channel_uniref50.txt_Kyte\ &\ Doolittle\ hydrophobiciity\ scale.png
    │   ├── 17_43_38__24_01_2018_voltage_gated_ion_channel_uniref50.txt_Eisenberg\ hydrophobiciity\ scale.pdf
    │   ├── 17_43_38__24_01_2018_voltage_gated_ion_channel_uniref50.txt_Eisenberg\ hydrophobiciity\ scale.png
    │   ├── 17_43_41__24_01_2018_voltage_gated_ion_channel_uniref50.txt_\ White\ and\ Wimley\ hydrophobiciity\ scale.pdf
    │   ├── 17_43_41__24_01_2018_voltage_gated_ion_channel_uniref50.txt_\ White\ and\ Wimley\ hydrophobiciity\ scale.png
    │   ├── 17_46_57__24_01_2018_voltage_gated_ion_channel_uniref50.txt_Complexity.pdf
    │   ├── 17_46_57__24_01_2018_voltage_gated_ion_channel_uniref50.txt_Complexity.png
    │   ├── 17_46_59__24_01_2018_voltage_gated_ion_channel_uniref50.txt_Length.pdf
    │   ├── 17_46_59__24_01_2018_voltage_gated_ion_channel_uniref50.txt_Length.png
    │   ├── 17_47_01__24_01_2018_voltage_gated_ion_channel_uniref50.txt_Kyte\ &\ Doolittle\ hydrophobiciity\ scale.pdf
    │   ├── 17_47_01__24_01_2018_voltage_gated_ion_channel_uniref50.txt_Kyte\ &\ Doolittle\ hydrophobiciity\ scale.png
    │   ├── 17_47_03__24_01_2018_voltage_gated_ion_channel_uniref50.txt_Eisenberg\ hydrophobiciity\ scale.pdf
    │   ├── 17_47_03__24_01_2018_voltage_gated_ion_channel_uniref50.txt_Eisenberg\ hydrophobiciity\ scale.png
    │   ├── 17_47_05__24_01_2018_voltage_gated_ion_channel_uniref50.txt_\ White\ and\ Wimley\ hydrophobiciity\ scale.pdf
    │   ├── 17_47_05__24_01_2018_voltage_gated_ion_channel_uniref50.txt_\ White\ and\ Wimley\ hydrophobiciity\ scale.png
    │   ├── 17_49_12__24_01_2018_voltage_gated_ion_channel_uniref50.txt_Complexity.pdf
    │   ├── 17_49_12__24_01_2018_voltage_gated_ion_channel_uniref50.txt_Complexity.png
    │   ├── 17_49_14__24_01_2018_voltage_gated_ion_channel_uniref50.txt_Length.pdf
    │   ├── 17_49_14__24_01_2018_voltage_gated_ion_channel_uniref50.txt_Length.png
    │   ├── 17_49_16__24_01_2018_voltage_gated_ion_channel_uniref50.txt_Kyte\ &\ Doolittle\ hydrophobiciity\ scale.pdf
    │   ├── 17_49_16__24_01_2018_voltage_gated_ion_channel_uniref50.txt_Kyte\ &\ Doolittle\ hydrophobiciity\ scale.png
    │   ├── 17_49_17__24_01_2018_voltage_gated_ion_channel_uniref50.txt_Eisenberg\ hydrophobiciity\ scale.pdf
    │   ├── 17_49_17__24_01_2018_voltage_gated_ion_channel_uniref50.txt_Eisenberg\ hydrophobiciity\ scale.png
    │   ├── 17_49_19__24_01_2018_voltage_gated_ion_channel_uniref50.txt_\ White\ and\ Wimley\ hydrophobiciity\ scale.pdf
    │   ├── 17_49_19__24_01_2018_voltage_gated_ion_channel_uniref50.txt_\ White\ and\ Wimley\ hydrophobiciity\ scale.png
    │   ├── 17_50_28__30_01_2018_KW\ porin.txt_Complexity.pdf
    │   ├── 17_50_28__30_01_2018_KW\ porin.txt_Complexity.png
    │   ├── 17_50_29__30_01_2018_KW\ porin.txt_Kyte\ &\ Doolittle\ hydrophobiciity\ scale.pdf
    │   ├── 17_50_29__30_01_2018_KW\ porin.txt_Kyte\ &\ Doolittle\ hydrophobiciity\ scale.png
    │   ├── 17_50_29__30_01_2018_KW\ porin.txt_Length.pdf
    │   ├── 17_50_29__30_01_2018_KW\ porin.txt_Length.png
    │   ├── 17_50_30__30_01_2018_KW\ porin.txt_Eisenberg\ hydrophobiciity\ scale.pdf
    │   ├── 17_50_30__30_01_2018_KW\ porin.txt_Eisenberg\ hydrophobiciity\ scale.png
    │   ├── 17_50_31__30_01_2018_KW\ porin.txt_\ White\ and\ Wimley\ hydrophobiciity\ scale.pdf
    │   ├── 17_50_31__30_01_2018_KW\ porin.txt_\ White\ and\ Wimley\ hydrophobiciity\ scale.png
    │   ├── 17_51_15__30_01_2018_KW\ porin.txt_Complexity.pdf
    │   ├── 17_51_15__30_01_2018_KW\ porin.txt_Complexity.png
    │   ├── 17_51_16__30_01_2018_KW\ porin.txt_Length.pdf
    │   ├── 17_51_16__30_01_2018_KW\ porin.txt_Length.png
    │   ├── 17_51_17__30_01_2018_KW\ porin.txt_Kyte\ &\ Doolittle\ hydrophobiciity\ scale.pdf
    │   ├── 17_51_17__30_01_2018_KW\ porin.txt_Kyte\ &\ Doolittle\ hydrophobiciity\ scale.png
    │   ├── 17_51_18__30_01_2018_KW\ porin.txt_\ White\ and\ Wimley\ hydrophobiciity\ scale.pdf
    │   ├── 17_51_18__30_01_2018_KW\ porin.txt_\ White\ and\ Wimley\ hydrophobiciity\ scale.png
    │   ├── 17_51_18__30_01_2018_KW\ porin.txt_Eisenberg\ hydrophobiciity\ scale.pdf
    │   ├── 17_51_18__30_01_2018_KW\ porin.txt_Eisenberg\ hydrophobiciity\ scale.png
    │   ├── 17_55_12__24_01_2018_KW\ ion\ channel.txt_Complexity.pdf
    │   ├── 17_55_12__24_01_2018_KW\ ion\ channel.txt_Complexity.png
    │   ├── 17_55_14__24_01_2018_KW\ ion\ channel.txt_Length.pdf
    │   ├── 17_55_14__24_01_2018_KW\ ion\ channel.txt_Length.png
    │   ├── 17_55_15__24_01_2018_KW\ ion\ channel.txt_Kyte\ &\ Doolittle\ hydrophobiciity\ scale.pdf
    │   ├── 17_55_15__24_01_2018_KW\ ion\ channel.txt_Kyte\ &\ Doolittle\ hydrophobiciity\ scale.png
    │   ├── 17_55_17__24_01_2018_KW\ ion\ channel.txt_Eisenberg\ hydrophobiciity\ scale.pdf
    │   ├── 17_55_17__24_01_2018_KW\ ion\ channel.txt_Eisenberg\ hydrophobiciity\ scale.png
    │   ├── 17_55_19__24_01_2018_KW\ ion\ channel.txt_\ White\ and\ Wimley\ hydrophobiciity\ scale.pdf
    │   ├── 17_55_19__24_01_2018_KW\ ion\ channel.txt_\ White\ and\ Wimley\ hydrophobiciity\ scale.png
    │   ├── 17_56_22__30_01_2018_KW\ porin.txt_Complexity.pdf
    │   ├── 17_56_22__30_01_2018_KW\ porin.txt_Complexity.png
    │   ├── 17_56_23__30_01_2018_KW\ porin.txt_Length.pdf
    │   ├── 17_56_23__30_01_2018_KW\ porin.txt_Length.png
    │   ├── 17_56_24__30_01_2018_KW\ porin.txt_Eisenberg\ hydrophobiciity\ scale.pdf
    │   ├── 17_56_24__30_01_2018_KW\ porin.txt_Eisenberg\ hydrophobiciity\ scale.png
    │   ├── 17_56_24__30_01_2018_KW\ porin.txt_Kyte\ &\ Doolittle\ hydrophobiciity\ scale.pdf
    │   ├── 17_56_24__30_01_2018_KW\ porin.txt_Kyte\ &\ Doolittle\ hydrophobiciity\ scale.png
    │   ├── 17_56_25__30_01_2018_KW\ porin.txt_\ White\ and\ Wimley\ hydrophobiciity\ scale.pdf
    │   ├── 17_56_25__30_01_2018_KW\ porin.txt_\ White\ and\ Wimley\ hydrophobiciity\ scale.png
    │   ├── 18_06_50__24_01_2018_KW\ ion\ channel.txt_Complexity.pdf
    │   ├── 18_06_50__24_01_2018_KW\ ion\ channel.txt_Complexity.png
    │   ├── 18_06_52__24_01_2018_KW\ ion\ channel.txt_Length.pdf
    │   ├── 18_06_52__24_01_2018_KW\ ion\ channel.txt_Length.png
    │   ├── 18_06_54__24_01_2018_KW\ ion\ channel.txt_Kyte\ &\ Doolittle\ hydrophobiciity\ scale.pdf
    │   ├── 18_06_54__24_01_2018_KW\ ion\ channel.txt_Kyte\ &\ Doolittle\ hydrophobiciity\ scale.png
    │   ├── 18_06_56__24_01_2018_KW\ ion\ channel.txt_Eisenberg\ hydrophobiciity\ scale.pdf
    │   ├── 18_06_56__24_01_2018_KW\ ion\ channel.txt_Eisenberg\ hydrophobiciity\ scale.png
    │   ├── 18_06_58__24_01_2018_KW\ ion\ channel.txt_\ White\ and\ Wimley\ hydrophobiciity\ scale.pdf
    │   ├── 18_06_58__24_01_2018_KW\ ion\ channel.txt_\ White\ and\ Wimley\ hydrophobiciity\ scale.png
    │   ├── 18_19_59__24_01_2018_KW\ ion\ channel.txt_Complexity.pdf
    │   ├── 18_19_59__24_01_2018_KW\ ion\ channel.txt_Complexity.png
    │   ├── 18_20_01__24_01_2018_KW\ ion\ channel.txt_Length.pdf
    │   ├── 18_20_01__24_01_2018_KW\ ion\ channel.txt_Length.png
    │   ├── 18_20_02__24_01_2018_KW\ ion\ channel.txt_Kyte\ &\ Doolittle\ hydrophobiciity\ scale.pdf
    │   ├── 18_20_02__24_01_2018_KW\ ion\ channel.txt_Kyte\ &\ Doolittle\ hydrophobiciity\ scale.png
    │   ├── 18_20_04__24_01_2018_KW\ ion\ channel.txt_Eisenberg\ hydrophobiciity\ scale.pdf
    │   ├── 18_20_04__24_01_2018_KW\ ion\ channel.txt_Eisenberg\ hydrophobiciity\ scale.png
    │   ├── 18_20_06__24_01_2018_KW\ ion\ channel.txt_\ White\ and\ Wimley\ hydrophobiciity\ scale.pdf
    │   ├── 18_20_06__24_01_2018_KW\ ion\ channel.txt_\ White\ and\ Wimley\ hydrophobiciity\ scale.png
    │   ├── 18_21_57__24_01_2018_KW\ ion\ channel.txt_Complexity.pdf
    │   ├── 18_21_57__24_01_2018_KW\ ion\ channel.txt_Complexity.png
    │   ├── 18_21_58__24_01_2018_KW\ ion\ channel.txt_Length.pdf
    │   ├── 18_21_58__24_01_2018_KW\ ion\ channel.txt_Length.png
    │   ├── 18_22_00__24_01_2018_KW\ ion\ channel.txt_Kyte\ &\ Doolittle\ hydrophobiciity\ scale.pdf
    │   ├── 18_22_00__24_01_2018_KW\ ion\ channel.txt_Kyte\ &\ Doolittle\ hydrophobiciity\ scale.png
    │   ├── 18_22_02__24_01_2018_KW\ ion\ channel.txt_Eisenberg\ hydrophobiciity\ scale.pdf
    │   ├── 18_22_02__24_01_2018_KW\ ion\ channel.txt_Eisenberg\ hydrophobiciity\ scale.png
    │   ├── 18_22_04__24_01_2018_KW\ ion\ channel.txt_\ White\ and\ Wimley\ hydrophobiciity\ scale.pdf
    │   ├── 18_22_04__24_01_2018_KW\ ion\ channel.txt_\ White\ and\ Wimley\ hydrophobiciity\ scale.png
    │   ├── 18_25_21__24_01_2018_KW\ ion\ channel.txt_Complexity.pdf
    │   ├── 18_25_21__24_01_2018_KW\ ion\ channel.txt_Complexity.png
    │   ├── 18_25_23__24_01_2018_KW\ ion\ channel.txt_Length.pdf
    │   ├── 18_25_23__24_01_2018_KW\ ion\ channel.txt_Length.png
    │   ├── 18_25_25__24_01_2018_KW\ ion\ channel.txt_Kyte\ &\ Doolittle\ hydrophobiciity\ scale.pdf
    │   ├── 18_25_25__24_01_2018_KW\ ion\ channel.txt_Kyte\ &\ Doolittle\ hydrophobiciity\ scale.png
    │   ├── 18_25_27__24_01_2018_KW\ ion\ channel.txt_Eisenberg\ hydrophobiciity\ scale.pdf
    │   ├── 18_25_27__24_01_2018_KW\ ion\ channel.txt_Eisenberg\ hydrophobiciity\ scale.png
    │   ├── 18_25_29__24_01_2018_KW\ ion\ channel.txt_\ White\ and\ Wimley\ hydrophobiciity\ scale.pdf
    │   ├── 18_25_29__24_01_2018_KW\ ion\ channel.txt_\ White\ and\ Wimley\ hydrophobiciity\ scale.png
    │   ├── 18_27_57__24_01_2018_KW\ ion\ channel.txt_Complexity.pdf
    │   ├── 18_27_57__24_01_2018_KW\ ion\ channel.txt_Complexity.png
    │   ├── 18_27_58__24_01_2018_KW\ ion\ channel.txt_Length.pdf
    │   ├── 18_27_58__24_01_2018_KW\ ion\ channel.txt_Length.png
    │   ├── 18_28_00__24_01_2018_KW\ ion\ channel.txt_Kyte\ &\ Doolittle\ hydrophobiciity\ scale.pdf
    │   ├── 18_28_00__24_01_2018_KW\ ion\ channel.txt_Kyte\ &\ Doolittle\ hydrophobiciity\ scale.png
    │   ├── 18_28_02__24_01_2018_KW\ ion\ channel.txt_Eisenberg\ hydrophobiciity\ scale.pdf
    │   ├── 18_28_02__24_01_2018_KW\ ion\ channel.txt_Eisenberg\ hydrophobiciity\ scale.png
    │   ├── 18_28_04__24_01_2018_KW\ ion\ channel.txt_\ White\ and\ Wimley\ hydrophobiciity\ scale.pdf
    │   ├── 18_28_04__24_01_2018_KW\ ion\ channel.txt_\ White\ and\ Wimley\ hydrophobiciity\ scale.png
    │   ├── 19_03_34__25_01_2018_7TM_list_UniRef50.txt_Complexity.pdf
    │   ├── 19_03_34__25_01_2018_7TM_list_UniRef50.txt_Complexity.png
    │   ├── 19_03_35__25_01_2018_7TM_list_UniRef50.txt_Length.pdf
    │   ├── 19_03_35__25_01_2018_7TM_list_UniRef50.txt_Length.png
    │   ├── 19_03_36__25_01_2018_7TM_list_UniRef50.txt_Eisenberg\ hydrophobiciity\ scale.pdf
    │   ├── 19_03_36__25_01_2018_7TM_list_UniRef50.txt_Eisenberg\ hydrophobiciity\ scale.png
    │   ├── 19_03_36__25_01_2018_7TM_list_UniRef50.txt_Kyte\ &\ Doolittle\ hydrophobiciity\ scale.pdf
    │   ├── 19_03_36__25_01_2018_7TM_list_UniRef50.txt_Kyte\ &\ Doolittle\ hydrophobiciity\ scale.png
    │   ├── 19_03_37__25_01_2018_7TM_list_UniRef50.txt_\ White\ and\ Wimley\ hydrophobiciity\ scale.pdf
    │   ├── 19_03_37__25_01_2018_7TM_list_UniRef50.txt_\ White\ and\ Wimley\ hydrophobiciity\ scale.png
    │   ├── 19_28_40__25_01_2018_7TM_list_UniRef50.txt_Complexity.pdf
    │   ├── 19_28_40__25_01_2018_7TM_list_UniRef50.txt_Complexity.png
    │   ├── 19_28_41__25_01_2018_7TM_list_UniRef50.txt_Length.pdf
    │   ├── 19_28_41__25_01_2018_7TM_list_UniRef50.txt_Length.png
    │   ├── 19_28_42__25_01_2018_7TM_list_UniRef50.txt_Kyte\ &\ Doolittle\ hydrophobiciity\ scale.pdf
    │   ├── 19_28_42__25_01_2018_7TM_list_UniRef50.txt_Kyte\ &\ Doolittle\ hydrophobiciity\ scale.png
    │   ├── 19_28_43__25_01_2018_7TM_list_UniRef50.txt_\ White\ and\ Wimley\ hydrophobiciity\ scale.pdf
    │   ├── 19_28_43__25_01_2018_7TM_list_UniRef50.txt_\ White\ and\ Wimley\ hydrophobiciity\ scale.png
    │   ├── 19_28_43__25_01_2018_7TM_list_UniRef50.txt_Eisenberg\ hydrophobiciity\ scale.pdf
    │   ├── 19_28_43__25_01_2018_7TM_list_UniRef50.txt_Eisenberg\ hydrophobiciity\ scale.png
    │   ├── 19_36_23__25_01_2018_7TM_list_UniRef50.txt_Complexity.pdf
    │   ├── 19_36_23__25_01_2018_7TM_list_UniRef50.txt_Complexity.png
    │   ├── 19_36_24__25_01_2018_7TM_list_UniRef50.txt_Length.pdf
    │   ├── 19_36_24__25_01_2018_7TM_list_UniRef50.txt_Length.png
    │   ├── 19_36_25__25_01_2018_7TM_list_UniRef50.txt_Kyte\ &\ Doolittle\ hydrophobiciity\ scale.pdf
    │   ├── 19_36_25__25_01_2018_7TM_list_UniRef50.txt_Kyte\ &\ Doolittle\ hydrophobiciity\ scale.png
    │   ├── 19_36_26__25_01_2018_7TM_list_UniRef50.txt_Eisenberg\ hydrophobiciity\ scale.pdf
    │   ├── 19_36_26__25_01_2018_7TM_list_UniRef50.txt_Eisenberg\ hydrophobiciity\ scale.png
    │   ├── 19_36_27__25_01_2018_7TM_list_UniRef50.txt_\ White\ and\ Wimley\ hydrophobiciity\ scale.pdf
    │   ├── 19_36_27__25_01_2018_7TM_list_UniRef50.txt_\ White\ and\ Wimley\ hydrophobiciity\ scale.png
    │   ├── 19_49_12__25_01_2018_7TM_list_UniRef50.txt_Complexity.pdf
    │   ├── 19_49_12__25_01_2018_7TM_list_UniRef50.txt_Complexity.png
    │   ├── 19_49_12__25_01_2018_7TM_list_UniRef50.txt_Length.pdf
    │   ├── 19_49_12__25_01_2018_7TM_list_UniRef50.txt_Length.png
    │   ├── 19_49_13__25_01_2018_7TM_list_UniRef50.txt_Kyte\ &\ Doolittle\ hydrophobiciity\ scale.pdf
    │   ├── 19_49_13__25_01_2018_7TM_list_UniRef50.txt_Kyte\ &\ Doolittle\ hydrophobiciity\ scale.png
    │   ├── 19_49_14__25_01_2018_7TM_list_UniRef50.txt_Eisenberg\ hydrophobiciity\ scale.pdf
    │   ├── 19_49_14__25_01_2018_7TM_list_UniRef50.txt_Eisenberg\ hydrophobiciity\ scale.png
    │   ├── 19_49_15__25_01_2018_7TM_list_UniRef50.txt_\ White\ and\ Wimley\ hydrophobiciity\ scale.pdf
    │   ├── 19_49_15__25_01_2018_7TM_list_UniRef50.txt_\ White\ and\ Wimley\ hydrophobiciity\ scale.png
    │   ├── 22_29_03__02_08_2018_KW\ ion\ channel.txt_TMSOC\ z-score.pdf
    │   ├── 22_29_03__02_08_2018_KW\ ion\ channel.txt_TMSOC\ z-score.png
    │   ├── 22_29_04__02_08_2018_KW\ ion\ channel.txt_Length.pdf
    │   ├── 22_29_04__02_08_2018_KW\ ion\ channel.txt_Length.png
    │   ├── 22_29_05__02_08_2018_KW\ ion\ channel.txt_Sequence\ Entropy.pdf
    │   ├── 22_29_05__02_08_2018_KW\ ion\ channel.txt_Sequence\ Entropy.png
    │   ├── 22_33_02__02_08_2018_KW\ ion\ channel.txt_TMSOC\ z-score.pdf
    │   ├── 22_33_02__02_08_2018_KW\ ion\ channel.txt_TMSOC\ z-score.png
    │   ├── 22_33_04__02_08_2018_KW\ ion\ channel.txt_Length.pdf
    │   ├── 22_33_04__02_08_2018_KW\ ion\ channel.txt_Length.png
    │   ├── 22_33_05__02_08_2018_KW\ ion\ channel.txt_Sequence\ Entropy.pdf
    │   ├── 22_33_05__02_08_2018_KW\ ion\ channel.txt_Sequence\ Entropy.png
    │   ├── 22_35_27__02_08_2018_KW\ ion\ channel.txt_TMSOC\ z-score.pdf
    │   ├── 22_35_27__02_08_2018_KW\ ion\ channel.txt_TMSOC\ z-score.png
    │   ├── 22_35_28__02_08_2018_KW\ ion\ channel.txt_Length.pdf
    │   ├── 22_35_28__02_08_2018_KW\ ion\ channel.txt_Length.png
    │   ├── 22_35_29__02_08_2018_KW\ ion\ channel.txt_Sequence\ Entropy.pdf
    │   ├── 22_35_29__02_08_2018_KW\ ion\ channel.txt_Sequence\ Entropy.png
    │   ├── 22_36_22__02_08_2018_KW\ ion\ channel.txt_TMSOC\ z-score.pdf
    │   ├── 22_36_22__02_08_2018_KW\ ion\ channel.txt_TMSOC\ z-score.png
    │   ├── 22_36_23__02_08_2018_KW\ ion\ channel.txt_Length.pdf
    │   ├── 22_36_23__02_08_2018_KW\ ion\ channel.txt_Length.png
    │   ├── 22_36_24__02_08_2018_KW\ ion\ channel.txt_Sequence\ Entropy.pdf
    │   ├── 22_36_24__02_08_2018_KW\ ion\ channel.txt_Sequence\ Entropy.png
    │   ├── 22_37_36__02_08_2018_KW\ ion\ channel.txt_TMSOC\ z-score.pdf
    │   ├── 22_37_36__02_08_2018_KW\ ion\ channel.txt_TMSOC\ z-score.png
    │   ├── 22_37_37__02_08_2018_KW\ ion\ channel.txt_Length.pdf
    │   ├── 22_37_37__02_08_2018_KW\ ion\ channel.txt_Length.png
    │   ├── 22_37_38__02_08_2018_KW\ ion\ channel.txt_Sequence\ Entropy.pdf
    │   └── 22_37_38__02_08_2018_KW\ ion\ channel.txt_Sequence\ Entropy.png
    ├── readme.md
    ├── scripts
    │   ├── TMSOC.pl
    │   ├── Uniprot_to_graphs_and_stats.py
    │   ├── computeSeqAAEntropyGrpIVL.pm
    │   ├── computeSeqAAHydrophobicity.pm
    │   ├── computeSummaryStatistics.pm
    │   ├── computeZscore.pm
    │   ├── generateTMclassification.pm
    │   └── getAAHydrophobicity.pm
    └── source_data_txt_files
        ├── 7TM_list_UniRef50.txt
        ├── 7tmrlist.txt
        ├── Ammonia\ transport\ KW-0924\ Uniref50.txt
        ├── GPCR_UniRef50.txt
        ├── GPCR_UniRef50_test.txt
        ├── GPCR_controlledvocabularyquery_UniRef50.txt
        ├── KW\ calcium\ channel.txt
        ├── KW\ chloride\ channel.txt
        ├── KW\ integrin.txt
        ├── KW\ ion\ channel.txt
        ├── KW\ ion\ transport.txt
        ├── KW\ ligand\ gated\ ion\ channel.txt
        ├── KW\ porin.txt
        ├── KW\ potassium\ channel.txt
        ├── KW\ protein\ transport.txt
        ├── KW\ sodium\ channel.txt
        ├── KW\ sugar\ transport.txt
        ├── KW\ voltage\ gated\ ion\ channel.txt
        ├── KW_GPCR_UniRef50.txt
        ├── Membrane\ attack\ complex\tKW-0473\ Uniref50.txt
        ├── Opsins_UniRef50.txt
        ├── Sulfate\ transport\tKW-0764\ Uniref50.txt
        ├── frizzledsmooth_UniRef50.txt
        ├── fungalmating_UniRef50.txt
        ├── ig_uniref50.txt
        ├── ion_channel_uniref50.txt
        ├── ligand_gated_ion_channel_uniref50.txt
        ├── metabotropicglutamate_UniRef50.txt
        ├── pore_UniRef50.txt
        ├── sugar_transporters_uniref50.txt
        ├── swissprot_transmem.txt
        ├── t2r_UniRef50.txt
        ├── tetraspanins.txt
        └── voltage_gated_ion_channel_uniref50.txt
