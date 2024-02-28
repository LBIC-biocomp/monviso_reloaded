# FileHandler Class

## Overview

The `FileHandler` class provides a simple interface for file and directory operations, including creation, deletion, reading, writing, and moving files. It's designed to unify the complexity of file handling in all the classes of the module, making it easier to fix bugs.

## Methods

### `remove_file(self, file_path: Union[str, Path])`

Removes a file at the specified path.

- **Parameters:**
  - `file_path`: The path to the file to be removed.

### `create_directory(self, dir_path: Union[str, Path])`

Creates a new directory at the specified path. If necessary, it will create all intermediate directories.

- **Parameters:**
  - `dir_path`: The path where the new directory will be created.

### `move_file(self, src: Union[str, Path], dest: Union[str, Path])`

Moves a file from the source path to the destination path.

- **Parameters:**
  - `src`: The source file path.
  - `dest`: The destination path.

### `copy_file(self, src: Union[str, Path], dest: Union[str, Path])`

Copies a file from the source path to the destination path.

- **Parameters:**
  - `src`: The source file path.
  - `dest`: The destination path.

### `write_file(self, file_path: Union[str, Path], content: str)`

Writes content to a file at the specified path. If the file exists, it will be removed before writing the new content.

- **Parameters:**
  - `file_path`: The file path where the content will be written.
  - `content`: The content to write to the file.

### `check_existence(self, path: Union[str, Path]) -> bool`

Checks whether a file or directory exists at the specified path.

- **Parameters:**
  - `path`: The path of the file or directory to check.

- **Returns:**
  - `bool`: `True` if the file or directory exists, `False` otherwise.

### `read_file(self, path: Union[str, Path]) -> str`

Reads and returns the content of a file at the specified path.

- **Parameters:**
  - `path`: The path of the file to read.

- **Returns:**
  - `str`: The content of the file.

