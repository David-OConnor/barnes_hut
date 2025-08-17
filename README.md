# Runs the Barnes Hut algorithm for fast n-body simulations.

[![Crate](https://img.shields.io/crates/v/barnes_hut.svg)](https://crates.io/crates/barnes_hut)
[![Docs](https://docs.rs/lin_alg/badge.svg)](https://docs.rs/barnes_hut)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15616833.svg)](https://doi.org/10.5281/zenodo.15616833)




![Visualization of cubes from this library](/tree.png)

This algorithm uses Tree Code to group source bodies, as an approximation. It leads to $O(n \log{} n)$ computation time, where $n$ is the number of bodies. Canonical use cases include gravity, and charged particle simulations.

This is much faster than a naive N-body approach $O(n^2)$, at high body counts. It is slower than Fast Multipole Methods (FMM), which group target bodies in addition to source ones. $O(n)$ Note that FMM are more difficult to implement; if you have a good implementation, use that instead of this library. There are currently, as of January 2025, no published implementations in Rust.

It uses the [Rayon library](https://docs.rs/rayon/latest/rayon/) to parallelize of the loop summing forces from each source on a given target. You may wish to also parallelize the loop over targets in your application code.

It uses the [lin_alg](https://crates.io/crates/lin_alg) library for the `Vec3` vector type. You will need to import this in your application code, and covert to it when implementing `BodyModel`.

Currently hard-coded for `f64`. Post an issue on GitHub if you'd like `f32` support.

This library is used to compute force or acceleration between pairs of bodies. It can be used to compute electric, or gravitational force,
for example. It can also be used to calculate gravitational acceleration directly, if set up as such using your `force` function.
this represents the gravitational mass of the target body cancelling with its inertial mass. $a=f/m$

If you're using periodic boundary conditions, for example to model solvents in structural biology,
this may not be appropriate; consider using SPME or similar. For example, in the [Ewald](https://github.com/david-oconnor/ewald)
lib.

It's generic over body type and force (or acceleration) function. Here's an example of adapting your type for use here:

```rust
use lin_alg::f64::Vec3;

impl BodyModel for Body {
    fn posit(&self) -> Vec3 {
        self.posit
    }

    fn mass(&self) -> f64 {
        // Or `self.charge`.
        self.mass
    }
}
```

Example use:

```rust
use barnes_hut::{self, BhConfig, Tree};
use rayon::prelude::*;

fn run_timesteps() {
    // Note: `BhConfig` Includes a `Default` implementation.
    let config = BhConfig {
        // The primary degree of freedom. 0 means no grouping. Higher values group more aggressively, leading to
        // less accurate, faster computation. 0.5 and 1.0 are common defaults.
        θ: 0.5,
        max_bodies_per_node: 1,
        // Safety, e.g. if bodies are very close together.
        max_tree_depth: 15,
    };
    
    for t in timesteps {
        // Create the tree of source bodies once per time step. Note that for this example, source and target
        // bodies are from the same set, but they don't have to be.
        let tree = Tree::new(&state.bodies, &bb, &state.config.bh_config);

        // Iterate, in parallel, over target bodies. The loop over source bodies is handled
        // by the acceleration function.
        bodies.par_iter_mut().enumerate().for_each(|(id, body_target)| {
            integrate(&config, &tree, target, id);
        });
    }
}

fn integrate(bh_config: &BhConfig, tree: &Tree, target: Body, id_target: usize) {
    // ...

    // This force or acceleration function can be whatever you'd like. This example shows Newtonian
    // Gravity with MOND.
    // `force_fn` accepts parameters (position vector, source mass or charge,  distance), and outputs
    // a force or acceleration vector.
    let force_fn = |acc_dir, mass_src, dist| {
        acc_newton_with_mond(acc_dir, mass, dist, Some(mond_fn), softening_factor_sq)
    };

    let accel = barnes_hut::run_bh(
        target.posit,
        id_target, // `id_target` is used to prevent self-interaction.
        tree,
        bh_config,
        &acc_fn,
    );
}
```