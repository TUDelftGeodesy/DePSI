"""Config class for depsi."""
from logging import Logger
import os
import yaml
from pydantic import BaseModel, ConfigDict
from pathlib import Path

logger = Logger(__name__)

class PSSelection(BaseModel):
    """
    Class for Point Selection parameters.
    """
    method: str
    threshold: float
    output_chunks: int
    model_config = ConfigDict(
        validate_default = True,
        populate_by_name = True,
        validate_assignment = True,
        extra = "forbid",
        )


class STMAddIncremental(BaseModel):
    """
    Class for STM Add Incremental parameters.
    """
    mode: str
    method: str
    recalibration_jump_size: float
    model_config = ConfigDict(extra = "forbid")


class STMPartitioning(BaseModel):
    """
    Class for STM Partitioning parameters.
    """
    db_partitioning: bool
    search_method: str
    cost_model: str
    min_partition_size: int
    model_config = ConfigDict(extra = "forbid")


class STMPartitioningVariables(BaseModel):
    """
    Class for STM Partitioning parameters.
    """
    amplitude_variable_name: str
    output_variable_prefix: str
    output_variables: list[str]
    model_config = ConfigDict(extra = "forbid")


class DetectOutliersSTM(BaseModel):
    """
    Class for Detect Outliers parameters.
    """
    db_outlier_detection: bool
    window_size: int
    n_sigma: int
    model_config = ConfigDict(extra = "forbid")


class Config(BaseModel):
    """
    Base class for configuration parameters.
    Common parameters and methods for all config classes.
    """
    name: str
    description: str
    paths: dict[str, str]
    model_config = ConfigDict(
        validate_default = True,
        populate_by_name = True,
        validate_assignment = True,
        extra = "allow",  # for testing, otherwise "forbid"
        )

    @classmethod
    def from_yaml(cls, config_file):
        """Read configs from a config.yaml file.

        If key is not found in config.yaml, the default value is used.
        """
        if not Path(config_file).exists():
            raise FileNotFoundError(f"Config file {config_file} not found.")

        with open(config_file, "r") as f:
            try:
                cfg = yaml.safe_load(f)
            except yaml.YAMLError as exc:
                raise SyntaxError(f"Error parsing config file {config_file}.") from exc
        return cls(**cfg)

    def to_yaml(self, config_file):
        """Write configs to a yaml config_file."""
        if Path(config_file).exists():
            logger.warning(f"Overwriting config file {config_file}.")

        cfg = self.model_dump(by_alias=True, warnings=False)
        with open(config_file, "w") as f:
            yaml.dump(cfg, f, sort_keys=False)


class ConfigGenerateSTM(Config):
    """
    Config class for generating STM.
    """
    name: str
    description: str
    ps_selection : PSSelection
    stm_add_incremental : STMAddIncremental
    project_stm_coordinates_projection: str
    stm_single_difference_mother: str
    do_stm_partitioning: bool
    stm_partitioning: STMPartitioning
    stm_partitioning_normal: STMPartitioningVariables
    stm_partitioning_single_difference: STMPartitioningVariables
    do_detect_outliers_stm: bool
    detect_outliers_stm: DetectOutliersSTM


def get_config_file():
    """Get the config file path."""
    config_path = Path(__file__).parent / ".config" / "depsi"
    if os.environ.get("CONFIG_PATH"):
        return os.environ.get("CONFIG_PATH")
    elif os.path.exists(config_path):
        yml_files = Path.glob(config_path, "*.yml")
        if len(yml_files) > 1:
            raise ValueError(
                f"Multiple config files found in {config_path}. Please specify one."
                )
        return config_path / yml_files[0]
    else:
        raise FileNotFoundError(
            "Config file not found."
            f"Please specify one in {config_path}"
            " or as an environment variable `CONFIG_PATH`."
        )