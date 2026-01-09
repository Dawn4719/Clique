## Accelerating $k$-Clique Listing via Vertex Classification and Repeated Branch Reduction
This repository provides the implementation of algorithms for $k$-clique listing problem.  
It includes our proposed methods VRkCL-V and VRkCL-E.

## Code Structure
- `VRkCL-V` – The source code of VRkCL-V, which is vertex-oriented $k$-clique listing algorithm.
- `VRkCL-E` – The source code of VRkCL-E, which is edge-oriented $k$-clique listing algorithm.
- `Dataset` – Example dataset.

---

## Compile the codes
When you already download the codes, run the following commands to compile our codes.
```
cd VRkCL-E
mkdir build
cd build
cmake ..
make
```
After running the codes, there will be executable files called `VRkCL-E` in `build/` directory.

---

## Datasets

- **Default datasets** (e.g., `WK`, `FB`) are included in `Dataset`.
- **Additional datasets** can be downloaded from [SNAP](https://snap.stanford.edu/data/index.html) and [Networkrepository](https://networkrepository.com/).
---

## Run the procedure
### Example Commands
- **Run 8-clique listing on wk dataset in `VRkCL-V`:**
```bash
./VRkCL-V wk 8
```
- **Run 8-clique listing on wk dataset in `VRkCL-E`:**
```bash
./VRkCL-E wk 8
```

---

## External Dependencies
- **BitCol** (baselines): available at [k-clique-listing](https://github.com/zer0y/k-clique-listing).
- **SDegree** (baselines): available at [k-clique-listing](https://github.com/zer0y/k-clique-listing).
- **CCR** (baselines): available at [coreCliqueRemoval](https://github.com/LightWant/coreCliqueRemoval).
- **EBBkC+ET** (baselines): available at [EBBkC](https://github.com/wangkaixin219/EBBkC).
---
**NOTE**: For a fair comparison, each algorithm must perform the full clique enumeration, combinatorial shortcuts are not permitted.



