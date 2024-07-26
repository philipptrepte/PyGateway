# PyGateway
PyGateway is a Python script to generate expression constructs by performing an in silico LR gateway reaction with entry and destination vectors. The script supports `.xdna`, `.dna` and `.gb` file formats and allows for the addition of features from a CSV file. The output vector files are in `.gb` file format.

**New entry vector** maps can be created using the [ccsb-vectormaps](https://github.com/csecker/ccsb-vectormaps) repository from [Christopher Secker](https://github.com/csecker).

## Prerequisites

- Python 3.x
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

### Example

```sh
python create_expression_vector.py --entry_vector_folder import/entry_vectors --dest_vector_folder import/destination_vectors --output_folder output --features_file import/features/all_features.csv
```

## Functions

- `read_xdna_file(xdna_file)`: Reads a `.xdna` file and transforms all sequences to upper case.

- `read_dna_file(dna_file)`: Reads a `.dna` file and transforms all sequences to upper case.

- `read_gb_file(gb_file)`: Reads a `.gb` file and transforms all sequences to upper case.

- `create_expression_vector(entry_vector_file, destination_vector_file, output_folder, features_file)`: Creates an expression vector from the given entry and destination vector files, and saves it to the specified output folder.

- `process_folder(entry_vector_folder, dest_vector_folder, features_file, output_folder)`: Parses the `create_expression_vector` function over the `entry_vector_folder` and `dest_vector_folder` folders to generate expression vectors with feature maps from the `feature_files` files and saves all expression vectors in the `output_folder`

## License

This project is licensed under the GPL-3.0 license - see the [LICENSE](LICENSE) file for details.

## Acknowledgments

- [Biopython](https://biopython.org/)
- [snapgene-reader](https://github.com/Edinburgh-Genome-Foundry/snapgene_reader)
- [GitHub Copilot](https://github.com/features/copilot)
