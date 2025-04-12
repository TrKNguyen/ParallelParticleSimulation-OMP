# Parallel Particle Collision Simulator - OpenMP

**1. Overview:**
This project implements a parallel particle collision simulator using OpenMP, designed to simulate the movement and collisions of particles in a 2D space. The simulation calculates particle interactions and updates velocities based on elastic collisions with walls and other particles.

**2. Algorithm:**
The algorithm uses a simple, brute-force approach to simulate the motion of particles. Each particle’s position is updated based on its velocity, and collisions are resolved based on predefined rules for particle-wall and particle-particle interactions. The simulation progresses in discrete timesteps, where velocities and positions are updated sequentially.

**3. Parallelization Strategy:**
- The simulation utilizes OpenMP for parallel execution.
- Particle interactions are computed in parallel, with collision detection handled within parallel loops.
- The program divides the task of updating particle positions and resolving collisions across multiple threads for improved performance.

**4. Optimization Efforts:**
- **First Optimization:** Parallelized the collision detection and resolution process, reducing runtime by distributing particle checks across multiple threads.
- **Second Optimization:** Implemented spatial partitioning to reduce the number of collision checks between distant particles, improving efficiency.

**5. Files Included:**
- **`sim.cc`:** Contains the core implementation of the particle simulation.
- **`sim_validator.h`:** Validation library for checking the correctness of the simulation.
- **`Makefile`:** Used to compile the OpenMP project.
- **`collision.h`:** Reference implementation for particle collision detection and resolution.
- **`Report.pdf`:** Detailed implementation report and performance analysis.
- **`viz.py`:** Python script to visualize the simulation as an `.mp4` file.
- **`run_bench.sh`:** Shell script for running the simulation on Slurm and benchmarking the performance.
- **`io.h`/`io.cc`:** Reference input/output code for reading and writing the simulation data.
- **`gen_testcase.py`:** Python script to generate random test cases for the simulation.
