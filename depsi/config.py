"""Config class for depsi."""

import os
from logging import Logger
from pathlib import Path

import yaml
from pydantic import BaseModel, ConfigDict
from pydantic_core import PydanticUndefined

logger = Logger(__name__)


class PSSelection(BaseModel):
    """Class for Point Selection parameters."""

    method: str
    threshold: float
    output_chunks: int
    model_config = ConfigDict(extra="forbid")


class STMAddIncremental(BaseModel):
    """Class for STM Add Incremental parameters."""

    mode: str
    method: str
    recalibration_jump_size: float
    model_config = ConfigDict(extra="forbid")


class STMPartitioning(BaseModel):
    """Class for STM Partitioning parameters."""

    db_partitioning: bool | None = False
    search_method: str
    cost_model: str
    min_partition_size: int
    model_config = ConfigDict(extra="forbid")


class STMPartitioningVariables(BaseModel):
    """Class for STM Partitioning parameters."""

    amplitude_variable_name: str
    output_variable_prefix: str
    output_variables: list[str]
    model_config = ConfigDict(extra="forbid")


class DetectOutliersSTM(BaseModel):
    """Class for Detect Outliers parameters."""

    db_outlier_detection: bool | None = True
    window_size: int
    n_sigma: int
    model_config = ConfigDict(extra="forbid")


class Config(BaseModel):
    """Base class for configuration parameters.

    Common parameters and methods for all config classes.
    """

    name: str = "basic_config"
    description: str = "Basic configuration."
    paths: dict[str, str] = {"working_dir": "."}
    model_config = ConfigDict(
        validate_default=True,
        populate_by_name=True,
        validate_assignment=True,
        extra="allow",
    )

    @classmethod
    def from_yaml(cls, config_file):
        """Read configs from a config.yaml file.

        If key is not found in config.yaml, the default value is used.
        """
        if not Path(config_file).exists():
            raise FileNotFoundError(f"Config file {config_file} not found.")

        with open(config_file) as f:
            try:
                cfg = yaml.safe_load(f)
            except yaml.YAMLError as exc:
                raise SyntaxError(f"Error parsing config file {config_file}.") from exc
        return cls(**cfg)

    @classmethod
    def to_yaml(cls, config_file):
        """Write configs to a yaml config_file."""
        if Path(config_file).exists():
            logger.warning(f"Overwriting config file {config_file}.")

        cfg = _schema(cls)
        with open(config_file, "w") as f:
            yaml.dump(cfg, f, sort_keys=False)


class ConfigGenerateSTM(Config):
    """Config class for generating STM."""

    name: str = "generate_stm"
    description: str = "Configuration for generating STM."
    ps_selection: PSSelection
    stm_add_incremental: STMAddIncremental
    project_stm_coordinates_projection: str
    stm_single_difference_mother: str
    do_stm_partitioning: bool = True
    stm_partitioning: STMPartitioning
    stm_partitioning_normal: STMPartitioningVariables
    stm_partitioning_single_difference: STMPartitioningVariables
    do_detect_outliers_stm: bool = True
    detect_outliers_stm: DetectOutliersSTM


def get_config_file():
    """Get the config file path."""
    config_path = Path.home() / ".config" / "depsi"
    if os.environ.get("CONFIG_PATH"):
        return os.environ.get("CONFIG_PATH")
    elif os.path.exists(config_path):
        yml_files = Path.glob(config_path, "*.yml")
        if len(yml_files) > 1:
            raise ValueError(f"Multiple config files found in {config_path}. Please specify one.")
        return config_path / yml_files[0]
    else:
        raise FileNotFoundError(
            f"Config file not found. Please specify one in {config_path} or as an environment variable `CONFIG_PATH`."
        )


def _schema(model: type[BaseModel]) -> dict:
    """Return the schema of a pydantic model."""
    schema_dict = {}
    for name, field in model.model_fields.items():
        field_type = field.annotation
        if isinstance(field_type, type) and issubclass(field_type, BaseModel):
            schema_dict[name] = _schema(field_type)
        else:
            schema_dict[name] = field.default
            if schema_dict[name] == PydanticUndefined:
                schema_dict[name] = f"undefined with type ({field_type.__name__})"
    return schema_dict
