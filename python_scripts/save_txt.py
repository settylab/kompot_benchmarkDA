import anndata
import scanpy as sc


def save_list_to_file(numbers, file_path, delimiter=''):
    """
    Saves a list of lists to a text file, with each inner list's elements separated by a specified delimiter.

    Parameters:
    - data: List of lists to be saved.
    - file_path: String specifying the path and name of the file to save the data to.
    - delimiter: String used to separate elements within each inner list. Defaults to ', '.
    """
    with open(file_path, 'w') as file:
        for number in numbers:
            # Write each number to the file, converting it to a string
            # and adding a newline character ('\n') after each number
            file.write(f"{number}\n")


def save_adata(adata,output_filename):
    adata.write(output_filename)
    