import os
import shutil
from typing import Any

class FileHandler:
    def __init__(self):
        pass

    def remove_file(self, file_path: str):
        """ Removes a file at the specified path. """
        try:
            os.remove(file_path)
            print(f"File {file_path} removed successfully.")
        except FileNotFoundError:
            print(f"The file {file_path} does not exist.")
        except Exception as e:
            print(f"Error removing file {file_path}: {e}")

    def create_directory(self, dir_path: str):
        """ Creates a new directory at the specified path. """
        try:
            os.makedirs(dir_path, exist_ok=True)
            print(f"Directory {dir_path} created successfully.")
        except Exception as e:
            print(f"Error creating directory {dir_path}: {e}")

    def move_file(self, src: str, dest: str):
        """ Moves a file from src to dest. """
        try:
            shutil.move(src, dest)
            print(f"File moved from {src} to {dest}.")
        except FileNotFoundError:
            print(f"The file {src} does not exist.")
        except Exception as e:
            print(f"Error moving file from {src} to {dest}: {e}")
