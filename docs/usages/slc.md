# SLC related methods

## Converting coregistered interferogram stack to SLC stack

After reading the coregistered interferogram stack and the mother SLC with `sarxarray`, the SLC stack can be reconstructed using the [`ifg_to_slc`](https://tudelftgeodesy.github.io/DePSI/api_reference/#depsi.slc.ifg_to_slc) method. The method takes the mother SLC and the coregistered interferogram stack as input and returns the reconstructed SLC stack.

```python
import sarxarray
from depsi.slc import ifg_to_slc
from pathlib import Path

f_mother_slc = 'path/to/mother_slc.raw' # Path to the mother SLC binary file
f_ifgs = list(sorted(Path('dir_ifgs').rglob("2*/ifnteferogram.raw")))  # List of paths of coregistered interferograms
shape = (10768, 40588) # Shape of the stack, (nrows, ncols)
reading_chunks = (2000, 2000)  # Reading chunks for lazy loading, (nrows, ncols)

# Lazy loading mother SLC and ifg stack
mother = sarxarray.from_binary([f_mother_slc], shape, dtype=np.complex64, chunks=reading_chunks)
ifgs = sarxarray.from_binary(f_ifgs, shape, dtype=np.complex64, chunks=reading_chunks)

# Generate reconstructed SLCs
slc_recon = ifg_to_slc(mother, ifgs)
```