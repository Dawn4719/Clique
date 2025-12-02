## Accelerating $k$-Clique Listing via Vertex Classification and Repeated Branch Reduction
This repository provides the implementation of algorithms for $k$-clique listing problem.  
It includes our proposed methods VR$k$CL-V and VR$k$CL-E.

## Code Structure

The main implementations are located in the `query` directory, with each algorithm in its own files:

- `BaseLine.hpp` – Baseline approach
- `IIT.hpp` – Interval Index Tree
- `IIT-R.hpp` – Roaring-compressed Interval Index Tree
- `DeltaGraph.h/cpp` – Delta-based temporal index
- `VILA.h/cpp` – Vectorized Interval Lifespan Algorithm
- `LLAMA.hpp` – Snapshot-based storage with page management
- `DELTA.h/cpp` – Variation-tolerant delta-based structure

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
- **Additional datasets** can be downloaded from [SNAP](https://snap.stanford.edu/data/index.html).
---

## Run the procedure
### Example Commands
- **Run $k$-clique listing `VRkCL-V`:**
```bash
./bin/VRkCL-V WK 8
```
- **Run $k$-clique listing `VRkCL-E`:**
```bash
./VRkCL-E WK 8
```

---

## External Dependencies
- **BitCol** (baselines): available at [k-clique-listing]https://github.com/zer0y/
- **SDegree** (baselines): available at [k-clique-listing]https://github.com/zer0y/k-clique-listing
- **CCR** (baselines): available at [coreCliqueRemoval](https://github.com/LightWant/coreCliqueRemoval)
- **EBBkC+ET** (baselines): available at [EBBkC](https://github.com/wangkaixin219/EBBkC)

---
