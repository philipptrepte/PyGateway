# PyGateway
PyGateway is a Python script to generate expression constructs by performing an in silico LR gateway reaction with entry and destination vectors. The script supports `.xdna`, `.dna` and `.gb` file formats and allows for the addition of features from a CSV file. The output vector files are in `.gb` file format and a unique identifier like `TEMP0001` will be preceding the file name, if not provided.

**New entry vector** maps can be created using the [ccsb-vectormaps](https://github.com/csecker/ccsb-vectormaps) repository from [Christopher Secker](https://github.com/csecker).

### Entry vector file naming

Entry vectors files need to start with `pDONR`, `pdonr`, `PDONR`, `pENTR`, `pentr`, `PENTR` followed by `201`, `221` or `223` followed by `_` and the `gene name`. 

An example is: `pDONR223_BAD_CCSB_3706.xdna`. Importantly, only `_` are accepted as name separators.

### Destination vector file naming

Destination vector files need to start with the backbone name, e.g. `pcDNA3.1` and `_` is the accepted name separator. Importantly, `GW` in the name will be replaced with the `gene name` of the entry vector. Additional information such as tags (`mCitrine`, `NL`) and tag site (`N1`, `C1`) can be kept if the respective information is separated in the file name with a `-`. 

For example: `pcDNA3.1_GW-mCitrine-N1.xdna` contains the information on the backbone `pcDNA3.1`, and the `GW` indicates the gene insertion site at the `N1`-terminus of the `mCitrine` protein.

## Prerequisites

- Python (was developed and tested using Python 3.12.4 on MacOS arm64)
- Biopython
- snapgene-reader

## Installation

1. Clone the repository:
    ```sh
    git clone https://github.com/philipptrepte/PyGateway.git
    cd PyGateway
    ```

2. Create conda environment and install required packages:
    ```sh
    conda create -n PyGateway
    conda activate PyGateway 
    conda config --add channels defaults
    conda config --add channels bioconda
    conda config --add channels conda-forge

    conda install biopython
    conda install snapgene-reader
    ```

## Usage

To create expression vector files, run the `create_expression_vector.py` script with the appropriate arguments:

```sh
python create_expression_vector.py --entry_vector_folder <entry_vector_folder> --dest_vector_folder <dest_vector_folder> --output_folder <output_folder> --features_file <features_file>
```

## Arguments

- `--entry_vector_folder` or `-e`: Path to the folder containing entry vector files.
- `--dest_vector_folder` or `-d`: Path to the folder containing destination vector files.
- `--output_folder` or `-o`: Path to the folder where the output files will be saved.
- `--features_file` or `-f`: Path to the CSV file containing features to be added.
- `--identifier` or `-i`: The identifier to use for the first expression vector. All following will be incremented by 1. Default to `TEMP0001`.

### Example

```sh
python create_expression_vector.py --entry_vector_folder import/entry_vectors --dest_vector_folder import/destination_vectors --output_folder output --features_file import/features/all_features.csv --identifier 
```

## Functions

- `read_xdna_file(xdna_file)`: Reads a `.xdna` file and transforms all sequences to upper case.

- `read_dna_file(dna_file)`: Reads a `.dna` file and transforms all sequences to upper case.

- `read_gb_file(gb_file)`: Reads a `.gb` file and transforms all sequences to upper case.

- `create_expression_vector(entry_vector_file, destination_vector_file, output_folder, features_file)`: Creates an expression vector from the given entry and destination vector files, and saves it to the specified output folder.

- `process_folder(entry_vector_folder, dest_vector_folder, features_file, output_folder)`: Parses the `create_expression_vector` function over the `entry_vector_folder` and `dest_vector_folder` folders to generate expression vectors with feature maps from the `feature_files` files and saves all expression vectors in the `output_folder`

## Conda Packages Information

| Name | Version | Build | Channel |
|------|---------|-------|---------|
| biopython | 1.84 | py312h7e5086c_0 | conda-forge |
| bzip2 | 1.0.8 | h99b78c6_7 | conda-forge |
| ca-certificates | 2024.7.4 | hf0a4a13_0 | conda-forge |
| html2text | 2024.2.26 | pyhd8ed1ab_1 | conda-forge |
| libblas | 3.9.0 | 22_osxarm64_openblas | conda-forge |
| libcblas | 3.9.0 | 22_osxarm64_openblas | conda-forge |
| libcxx | 18.1.8 | h167917d_0 | conda-forge |
| libexpat | 2.6.2 | hebf3989_0 | conda-forge |
| libffi | 3.4.2 | h3422bc3_5 | conda-forge |
| libgfortran | 5.0.0 | 13_2_0_hd922786_3 | conda-forge |
| libgfortran5 | 13.2.0 | hf226fd6_3 | conda-forge |
| liblapack | 3.9.0 | 22_osxarm64_openblas | conda-forge |
| libopenblas | 0.3.27 | openmp_h517c56d_1 | conda-forge |
| libsqlite | 3.46.0 | hfb93653_0 | conda-forge |
| libzlib | 1.3.1 | hfb2fe0b_1 | conda-forge |
| llvm-openmp | 18.1.8 | hde57baf_0 | conda-forge |
| ncurses | 6.5 | hb89a1cb_0 | conda-forge |
| numpy | 2.0.0 | py312hb544834_0 | conda-forge |
| openssl | 3.3.1 | hfb2fe0b_2 | conda-forge |
| pip | 24.0 | pyhd8ed1ab_0 | conda-forge |
| python | 3.12.4 | h30c5eda_0_cpython | conda-forge |
| python_abi | 3.12 | 4_cp312 | conda-forge |
| readline | 8.2 | h92ec313_1 | conda-forge |
| setuptools | 71.0.4 | pyhd8ed1ab_0 | conda-forge |
| snapgene-reader | 0.1.20 | pyh5e36f6f_0 | bioconda |
| tk | 8.6.13 | h5083fa2_1 | conda-forge |
| tzdata | 2024a | h0c530f3_0 | conda-forge |
| wheel | 0.43.0 | pyhd8ed1ab_1 | conda-forge |
| xmltodict | 0.13.0 | pyhd8ed1ab_0 | conda-forge |
| xz | 5.2.6 | h57fd34a_0 | conda-forge |


## License

This project is licensed under the GPL-3.0 license - see the [LICENSE](LICENSE) file for details.

## Acknowledgments

- [Biopython](https://biopython.org/)
- [snapgene-reader](https://github.com/Edinburgh-Genome-Foundry/snapgene_reader)
- [GitHub Copilot](https://github.com/features/copilot)
