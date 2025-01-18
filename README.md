# Runs the Barnes Hut algorithm for fast n-body simulations.

[![Crate](https://img.shields.io/crates/v/barnes_hut.svg)](https://crates.io/crates/barnes_hut)
[![Docs](https://docs.rs/lin_alg/badge.svg)](https://docs.rs/barnes_hut)

This algorithm uses Tree Code to group source bodies, as an approximation. It leads to O(N(log N)) computation time, where `N` is the number of bodies. Canonical use cases include gravity, and charged particle simulations.

This is much faster than a naive N-body approach, at high body counts. It is slower than Fast Multiple Methods (FMM),
which group target bodies in addition to source ones. Note that FMM methods are more difficult to implement, but if you have a good implementation, use that instead of this library. There are currently, as of January 2025, no published implementations in Rust.

It uses the [Rayon library](https://docs.rs/rayon/latest/rayon/) for parallelization of the loop summing forces from each source on a given target. You may wish to parallelize the loop over targets in your application code.

It usses the [lin_alg](https://crates.io/crates/lin_alg) library for the `Vec3` vector type. You will need to import this in your application code, and covert to it when implementing `BodyModel`.

Currently hard-coded for `f64`. Post an issue on Github if you'd like `f32` support.

This library is generic over `Body` and acceleration function . Here's an example of adapting your type for use here:

```rust
use lin_alg::f64::Vec3;

impl BodyModel for Body {
    fn posit(&self) -> Vec3 {
        // You may need to convert to the `lin_alg::f64::Vec3` type here, e.g:
        // Vec3::new(self.x, self.y, self.z)
        self.posit
    }

    fn mass(&self) -> f64 {
        self.mass
    }
}
```

Example use:

```rust
use barnes_hut::{self, BhConfig, Tree};

fn run_timestep() {
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
        // Create the tree once per time step.
        let tree = Tree::new(&state.bodies, &bb, &state.config.bh_config);

        // ...

        // Iterate, in parallel, over target bodies. The loop over source bodies, per target, is handled
        // by the acceleration function.
        // Or parallel:
        // bodies.par_iter_mut().enumerate().for_each(|(id, body_target)| { // ...
        for (id, target) in bodies.iter_mut().enumerate() {
            integrate(&config, &tree, target);
        }
    }
}

fn integrate(bh_config: &BhConfig, tree: &Tree, target: Body) {
    // ...

    // This acceleration function can be whatever you'd like. This example shows Newtonian
    // Gravity with MOND.
    let acc_fn = |acc_dir, mass, dist| {
        acc_newton_with_mond(acc_dir, mass, dist, Some(mond_fn), softening_factor_sq)
    };

    let accel = barnes_hut::acc_bh(
        target.posit,
        id_target, // Currently unused; set to 0 if you'd like, or the enumeration of bodies.
        tree,
        bh_config,
        &acc_fn,
    );
}
```