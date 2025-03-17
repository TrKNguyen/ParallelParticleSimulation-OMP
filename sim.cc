#include <omp.h>
//#include<iostream>
//#include<ctime>
#include <cmath>
#include <fstream>
#include <vector>
#include <utility>

#include "collision.h"
#include "io.h"
#include "sim_validator.h"

inline bool resolve_wall_collisions(std::vector<std::vector<Particle*>>& grid, int num_grid, int square_size, int radius) {
    bool collisionOccurred = false;
    #pragma omp for collapse(2)
    for (int r = 0; r < num_grid; r++) {
        for (int c = 0; c < num_grid; c++) {
            if (r != 0 && c != 0 && c != num_grid - 1 && r != num_grid - 1) continue;
            const std::vector<Particle*>& cell = grid[r * num_grid + c];
            for (Particle* p : cell) {
                if (is_wall_collision(p->loc, p->vel, square_size, radius)) {
                    resolve_wall_collision(p->loc, p->vel, square_size, radius);
                    collisionOccurred = true;
                }
            }      
        }
    }
    return collisionOccurred;
}

inline bool resolve_adjacent_collisions(std::vector<std::pair<Particle*, Particle*>>& candidate) {
    bool collisionOccurred = false;
    
    #pragma omp parallel for
    for (size_t i = 0; i < candidate.size(); ++i) {
        Particle* p1 = candidate[i].first;
        Particle* p2 = candidate[i].second;
        if (is_particle_moving_closer(p1->loc, p1->vel, p2->loc, p2->vel)) {
            resolve_particle_collision(p1->loc, p1->vel, p2->loc, p2->vel);
            collisionOccurred = true;
        }
    }
    return collisionOccurred;
}

void build_possible_collisions(std::vector<std::vector<Particle*>>& grid, std::vector<std::pair<Particle*, Particle*>>& candidate, int num_grid, int radius, int start_x, int start_y) {
    #pragma omp parallel for collapse(2) 
    for (int r = start_x; r < num_grid; r += 2) {
        for (int c = start_y; c < num_grid; c += 2) {
            int index = r * num_grid + c;
            std::vector<Particle*>& cell = grid[index];
            for (int dr = -1; dr <= 1; dr++) {
                for (int dc = 0; dc <= 1; dc++) {
                    if (dr == -1 && dc == 0) continue;
                    int nr = r + dr;
                    int nc = c + dc;
                    if (nr >= 0 && nr < num_grid && nc >= 0 && nc < num_grid) {
                        int neighborIndex = nr * num_grid + nc;
                        std::vector<Particle*>& neighborCell = grid[neighborIndex];
                        for (Particle* p1 : cell) {
                            for (Particle* p2 : neighborCell) {
                                if (is_particle_overlap(p1->loc, p2->loc, radius)) {
                                    #pragma omp critical
                                    candidate.emplace_back(p1, p2);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

void resolve_collisions(std::vector<std::vector<Particle*>>& grid, int num_grid, int square_size, int radius) {
    bool collisionOccurred = true;
    std::vector<std::pair<Particle*, Particle*>> candidate[4];

    #pragma omp parallel for collapse(2)
    for (int i = 0; i <= 1; i++) {
        for (int j = 0; j <= 1; j++) {
            build_possible_collisions(grid, candidate[(i << 1) | j], num_grid, radius, i, j);
        }
    }
    while (collisionOccurred) {
        collisionOccurred = false;
        collisionOccurred |= resolve_wall_collisions(grid, num_grid, square_size, radius);
        
        collisionOccurred |= resolve_adjacent_collisions(candidate[0]);
        collisionOccurred |= resolve_adjacent_collisions(candidate[1]);
        collisionOccurred |= resolve_adjacent_collisions(candidate[2]);
        collisionOccurred |= resolve_adjacent_collisions(candidate[3]);
    }
}

void simulate_step(std::vector<Particle>& particles, const Params& params) {
    int grid_size = 4 * params.param_radius;
    int num_grid = params.square_size / grid_size;
    std::vector<std::vector<Particle*>> grid(num_grid * num_grid);
    //std::cerr << "Number of grid each: " << num_grid << " Grid size: " << grid_size << std::endl;
    for (size_t i = 0; i < particles.size(); i++) {
        particles[i].loc.x += particles[i].vel.x;
        particles[i].loc.y += particles[i].vel.y;

        int idr = particles[i].loc.x / grid_size;
        int idc = particles[i].loc.y / grid_size;
        if (idr < 0) idr = 0;
        if (idr >= num_grid) idr = num_grid - 1;
        if (idc < 0) idc = 0;
        if (idc >= num_grid) idc = num_grid - 1;
        int index = idr * num_grid + idc;
        grid[index].push_back(&particles[i]);
    }
    resolve_collisions(grid, num_grid, params.square_size, params.param_radius);
}

int main(int argc, char* argv[]) {
    // Read arguments and input file
    Params params{};
    std::vector<Particle> particles;
    read_args(argc, argv, params, particles);

    // Set number of threads
    omp_set_num_threads(params.param_threads);

#if CHECK == 1
    // Initialize collision checker
    SimulationValidator validator(params.param_particles, params.square_size, params.param_radius, params.param_steps);
    // Initialize with starting positions
    validator.initialize(particles);
    // Uncomment the line below to enable visualization (makes program much slower)
    // validator.enable_viz_output("test.out");
#endif

    // TODO: this is the part where you simulate particle behavior.

    /*
    After simulating each timestep, you must call this exact block below.
    Make sure that your final submission has both the validation logic above and below included, within the #if

    #if CHECK == 1
        validator.validate_step(particles);
    #endif
    */
    for (int step = 0; step < params.param_steps; ++step) {
        simulate_step(particles, params);
        //std::cerr << "Current step: " << step << " " << (clock() - start) / CLOCKS_PER_SEC << std::endl;
        #if CHECK == 1
            validator.validate_step(particles);
        #endif
    }
}