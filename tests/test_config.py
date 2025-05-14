"""Test configuration loading and validation."""

import pytest
import yaml

from depsi.config import Config, ConfigGenerateSTM, get_config_file


@pytest.fixture
def stm_config_params():
    return {
        "name": "test_config",
        "description": "Test configuration for STM generation.",
        "paths": {
            "slc_file": "./test_slc.zarr",
            "output_dir": "./output",
        },
        "ps_selection": {
            "method": "test_method",
            "threshold": 0.5,
            "output_chunks": 100,
        },
        "stm_add_incremental": {
            "mode": "test_mode",
            "method": "test_method",
            "recalibration_jump_size": 10,
        },
        "project_stm_coordinates_projection": "RD",
        "stm_single_difference_mother": "auto",
        "do_stm_partitioning": True,
        "stm_partitioning": {
            "db_partitioning": False,
            "search_method": "test_search",
            "cost_model": "test_cost",
            "min_partition_size": 10,
        },
        "stm_partitioning_normal": {
            "amplitude_variable_name": "test_amplitude",
            "output_variable_prefix": "test_prefix",
            "output_variables": ["var1", "var2"],
        },
        "stm_partitioning_single_difference": {
            "amplitude_variable_name": "test_amplitude_sd",
            "output_variable_prefix": "test_prefix_sd",
            "output_variables": ["var1_sd", "var2_sd"],
        },
        "do_detect_outliers_stm": True,
        "detect_outliers_stm": {
            "db_outlier_detection": True,
            "window_size": 5,
            "n_sigma": 3,
        },
    }


@pytest.fixture
def basic_config_params():
    """Create a temporary YAML config file."""
    return {
        "name": "test_config",
        "description": "Test configuration.",
        "paths": {
            "work_dir": ".",
        },
    }


def write_config_file(params, dir):
    """Create a temporary YAML config file."""
    config_path = dir / "test_config.yaml"
    with open(config_path, "w") as f:
        yaml.dump(params, f, sort_keys=False)
    return config_path


class TestConfig:
    def test_config_loading(self, basic_config_params, tmp_path):
        """Test loading configuration from a YAML file."""
        cofig_path = write_config_file(basic_config_params, tmp_path)
        config = Config.from_yaml(cofig_path)
        assert config.name == "test_config"
        assert config.description == "Test configuration."
        assert config.paths["work_dir"] == "."

    def test_invalid_config_file(self):
        """Test loading configuration from a non-existent file."""
        with pytest.raises(FileNotFoundError):
            Config.from_yaml("invalid_config.yaml")

    def test_to_yaml(self, tmp_path):
        """Test saving configuration template to a YAML file."""
        config_path = tmp_path / "config_template.yml"
        Config.to_yaml(config_path)
        with open(config_path) as f:
            loaded_cfg = yaml.safe_load(f)
        assert loaded_cfg["name"] == "basic_config"
        assert loaded_cfg["description"] == "Basic configuration."
        assert loaded_cfg["paths"]["working_dir"] == "."


class TestSTMConfig:
    def test_config_loading(self, stm_config_params, tmp_path):
        """Test loading configuration from a YAML file."""
        config_path = write_config_file(stm_config_params, tmp_path)
        config = ConfigGenerateSTM.from_yaml(config_path)
        assert config.name == "test_config"
        assert config.description == "Test configuration for STM generation."
        assert config.paths["slc_file"] == "./test_slc.zarr"
        assert config.ps_selection.method == "test_method"
        assert config.do_stm_partitioning is True

    def test_to_yaml(self, tmp_path):
        """Test saving configuration template to a YAML file."""
        config_path = tmp_path / "config_template.yml"
        ConfigGenerateSTM.to_yaml(config_path)
        with open(config_path) as f:
            loaded_cfg = yaml.safe_load(f)
        assert loaded_cfg["name"] == "generate_stm"
        assert loaded_cfg["description"] == "Configuration for generating STM."
        assert loaded_cfg["paths"]["working_dir"] == "."
        assert loaded_cfg["ps_selection"]["method"] == "undefined with type (str)"

    def test_invalid_config(self, stm_config_params, tmp_path):
        """Test loading configuration with missing required fields."""
        invalid_config = stm_config_params.copy()
        del invalid_config["stm_single_difference_mother"]  # missing required field
        config_path = write_config_file(invalid_config, tmp_path)
        with pytest.raises(ValueError):
            ConfigGenerateSTM.from_yaml(config_path)

    def test_config_with_extra_allow(self, stm_config_params, tmp_path):
        """Test loading configuration with extra fields."""
        invalid_config = stm_config_params.copy()
        invalid_config["extra_field"] = "extra_value"  # Extra field
        config_path = write_config_file(invalid_config, tmp_path)
        config = ConfigGenerateSTM.from_yaml(config_path)
        assert config.extra_field == "extra_value"

    def test_config_with_invalid_field_type(self, stm_config_params, tmp_path):
        """Test loading configuration with invalid field type."""
        invalid_config = stm_config_params.copy()
        invalid_config["ps_selection"]["method"] = 123  # Invalid type
        config_path = write_config_file(invalid_config, tmp_path)
        with pytest.raises(ValueError):
            ConfigGenerateSTM.from_yaml(config_path)

    def test_config_with_extra_forbid(self, stm_config_params, tmp_path):
        """Test loading configuration with forbidden key in subclasses."""
        invalid_config = stm_config_params.copy()
        invalid_config["ps_selection"]["extra"] = 123  # forbid extra
        config_path = write_config_file(invalid_config, tmp_path)
        with pytest.raises(ValueError):
            ConfigGenerateSTM.from_yaml(config_path)


def test_get_config_file_env_var():
    """Test the get_config_file function."""
    # add a dummy path to env variable
    import os

    os.environ["CONFIG_PATH"] = "/dummy/path/to/config"

    config_path = get_config_file()
    assert config_path == "/dummy/path/to/config"


def test_get_config_file_no_default():
    """Test the get_config_file function."""
    # remove the env variable
    import os

    if "CONFIG_PATH" in os.environ:
        del os.environ["CONFIG_PATH"]

    with pytest.raises(FileNotFoundError):
        get_config_file()
