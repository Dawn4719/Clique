## Accelerating $k$-Clique Listing via Vertex Classification and Repeated Branch Reduction
This repository provides the implementation of algorithms for $k$-clique listing problem.  
It includes our proposed methods VR$k$CL-V and VR$k$CL-E.

## Code Structure
- `VRkCL-V` – The sorse code of VR$k$CL-V
- `VRkCL-E` – The sorse code of VR$k$CL-E
- `Dataset` – Example dataset

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
- **Run $k$-clique listing `VRkCL-V`:**
```bash
./VRkCL-V WK 8
```
- **Run $k$-clique listing `VRkCL-E`:**
```bash
./VRkCL-E WK 8
```

---

## External Dependencies
- **BitCol** (baselines): available at [k-clique-listing](https://github.com/zer0y/k-clique-listing)
- **SDegree** (baselines): available at [k-clique-listing](https://github.com/zer0y/k-clique-listing)
- **CCR** (baselines): available at [coreCliqueRemoval](https://github.com/LightWant/coreCliqueRemoval)
- **EBBkC+ET** (baselines): available at [EBBkC](https://github.com/wangkaixin219/EBBkC)

---
