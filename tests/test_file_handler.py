import os
from pathlib import Path
import pytest
from monviso_reloaded.file_handler import FileHandler

@pytest.fixture
def handler():
    with FileHandler() as fh:
        yield fh

def test_create_directory(handler, tmp_path):
    new_dir = tmp_path / "new_folder"
    assert not new_dir.exists()
    handler.create_directory(new_dir)
    assert new_dir.exists()
    assert new_dir.is_dir()

def test_remove_file(handler, tmp_path):
    f = tmp_path / "temp.txt"
    f.write_text("to be deleted")
    assert f.exists()
    handler.remove_file(f)
    assert not f.exists()

def test_write_and_read_file(handler, tmp_path):
    f = tmp_path / "file.txt"
    content = "hello world"
    handler.write_file(f, content)
    read_content = handler.read_file(f)
    assert read_content == content

def test_check_existence(handler, tmp_path):
    existing_file = tmp_path / "existing.txt"
    non_existing_file = tmp_path / "nope.txt"
    existing_file.write_text("hi")
    assert handler.check_existence(existing_file)
    assert not handler.check_existence(non_existing_file)

def test_move_file(handler, tmp_path):
    src = tmp_path / "source.txt"
    dest = tmp_path / "moved.txt"
    src.write_text("data")
    handler.move_file(src, dest)
    assert dest.exists()
    assert not src.exists()
    assert dest.read_text() == "data"

def test_copy_file(handler, tmp_path):
    src = tmp_path / "source.txt"
    dest = tmp_path / "copy.txt"
    src.write_text("copy this")
    handler.copy_file(src, dest)
    assert dest.exists()
    assert src.exists()
    assert dest.read_text() == "copy this"

def test_get_date(handler, tmp_path):
    f = tmp_path / "file.txt"
    f.write_text("check time")
    ctime = handler.get_date(f)
    assert isinstance(ctime, float)

def test_rename_files_in_directory(handler, tmp_path):
    f1 = tmp_path / "old_1.txt"
    f2 = tmp_path / "old_2.txt"
    f1.write_text("a")
    f2.write_text("b")
    
    handler.rename_files_in_directory(tmp_path, "old", "new")
    
    assert not (tmp_path / "old_1.txt").exists()
    assert not (tmp_path / "old_2.txt").exists()
    assert (tmp_path / "new_1.txt").exists()
    assert (tmp_path / "new_2.txt").exists()
