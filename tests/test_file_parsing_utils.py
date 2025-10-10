import pytest
from FetoFlow.file_parsing_utils import read_nodes, read_elements, define_fields_from_files


def test_read_nodes_invalid_type():
    with pytest.raises(TypeError):
        read_nodes(123)


def test_read_nodes_invalid_extension(tmp_path):
    f = tmp_path / "test.txt"
    f.write_text("not an ipnode file")
    with pytest.raises(TypeError):
        read_nodes(str(f))


def test_read_elements_invalid_type():
    with pytest.raises(TypeError):
        read_elements(123)


def test_read_elements_invalid_extension(tmp_path):
    f = tmp_path / "test.txt"
    f.write_text("not an ipelem file")
    with pytest.raises(TypeError):
        read_elements(str(f))


def test_define_fields_from_files_type_errors():
    with pytest.raises(TypeError):
        define_fields_from_files(123)
