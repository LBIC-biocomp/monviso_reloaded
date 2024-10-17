import os
import shutil
from pathlib import Path
from typing import Union

class FileHandler:
    def __init__(self):
        pass

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        pass

    def remove_file(self, file_path: Union[str, Path], verbose=False):
        """Removes a file at the specified path."""
        file_path = str(file_path)
        try:
            os.remove(file_path)
            if verbose:
                print(f"File {file_path} removed successfully.")
        except FileNotFoundError:
                print(f"The file {file_path} does not exist.")
        except Exception as e:
            print(f"Error removing file {file_path}: {e}")

    def create_directory(self, dir_path: Union[str, Path], verbose=False):
        """Creates a new directory at the specified path.
        It creates all the directory tree if it does not exist.
        """
        try:
            os.makedirs(str(dir_path), exist_ok=True)
            if verbose:
                print(f"Directory {str(dir_path)} created successfully.")
        except Exception as e:
            print(f"Error creating directory {dir_path}: {e}")

    def move_file(self, src: Union[str, Path], dest: Union[str, Path], verbose=False):
        """Moves a file from src to dest."""
        src = str(src)
        dest = str(dest)
        try:
            shutil.move(src, dest)
            if verbose:
                print(f"File moved from {src} to {dest}.")
        except FileNotFoundError:
            print(f"The file {src} does not exist.")
        except Exception as e:
            print(f"Error moving file from {src} to {dest}: {e}")

    def copy_file(self, src: Union[str, Path], dest: Union[str, Path],verbose=False):
        """Copies a file from src to dest."""
        src = str(src)
        dest = str(dest)
        try:
            shutil.copy(src, dest)
            if verbose:
                print(f"File copied from {src} to {dest}.")
        except FileNotFoundError:
            print(f"The file {src} does not exist.")
        except Exception as e:
            print(f"Error copying file from {src} to {dest}: {e}")

    def write_file(self, file_path: Union[str, Path], content: str, verbose=False):
        """Writes content to a file, removing it first if it exists."""

        # Check if the file exists and remove it
        # using the class's remove_file method
        file_path = str(file_path)
        if os.path.exists(file_path):
            self.remove_file(file_path,verbose)

        # Now, write the new content to the file
        try:
            with open(file_path, "w") as file:
                file.write(content)
            if verbose:
                print(f"Content written to {file_path} successfully.")
        except Exception as e:
            print(f"Error writing to file {file_path}: {e}")

    def check_existence(self, path: Union[str, Path]) -> bool:
        """
        Checks whether a given file or directory exists.

        Args:
            path (Union[str, Path]): The path of the file or directory
            to check.

        Returns:
            bool: True if the file or directory exists, False otherwise.
        """
        # Convert Path to string if it's not already a string
        path_str = str(path)
        return os.path.exists(path_str)

    def read_file(self, path: Union[str, Path]) -> str:
        """
        Read and return the content of a file after checking its existence.

        Args:
            path (Union[str, Path]): The path of the file or directory
            to check.

        Returns:
            str: The file content.
        """
        if not self.check_existence(path):
            raise (FileNotFoundError("The specified file was not found."))

        else:
            path = Path(path)
            with path.open("r") as file:
                content = file.read()
                return content
        
    def get_date(self,path: Union[str,Path]):
        """ Given a file, returns the creation a last-modified
            time as a tuple in "seconds since the epoch"

        Args:
            path (Union[str,Path]): The path of the file
        """
        # Both the variables would contain time
        # elapsed since EPOCH in float
        ti_c = os.path.getctime(path)
        ti_m = os.path.getmtime(path)

        return ti_c