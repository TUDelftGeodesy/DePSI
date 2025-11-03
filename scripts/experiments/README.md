# Configuration file in Depsi

In this folder, the python script `generate_stm_nl_amsterdam_s1_dsc_t037.py` implements several fucntions to generate a STM file. The script reads the input arguments of the functions from a yaml file `generate_stm_nl_amsterdam_s1_dsc_t037.yml` i.e. the configuration file.

## How to generate a config template in yml format

To generate a config template in yml format, in python:

```python
from depsi.config import ConfigGenerateSTM

ConfigGenerateSTM.to_yaml("your_config_file.yml")
```

The `your_config_file.yml` is the name of the config file you want to generate.
The generated config file contains all the parameters needed to generate a STM
file. You can modify the values of the parameters in the config file according
to your needs.

## How to read the config file

To read (and validate) the config file, in python:

```python
from depsi.config import ConfigGenerateSTM

config = ConfigGenerateSTM.from_yaml("your_config_file.yml")
print(config)
```

The `config` object contains all the parameters needed to run the functions,
generating a STM file. You can access a specific parameter by its key. For
example, to access the method parameter for `ps_selection`, you can do:

```python
print(config.ps_selection.method)
```

## How the config parameters are structured whithin template

There is a `config.py` module in Depsi that contains the `ConfigGenerateSTM`
class. This class is used to generate a config template in yml format and
also to read (and validate) parameters from a config file in yml format.
