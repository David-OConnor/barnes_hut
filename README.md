# Runs the Barnes Hut algorithm for fast n-body simulations.

[![Crate](https://img.shields.io/crates/v/barnes_hut.svg)](https://crates.io/crates/barnes_hut)
[![Docs](https://docs.rs/lin_alg/badge.svg)](https://docs.rs/barnes_hut)

This algorithm uses Tree Code to group source bodies, as an approximation. It leads to $O(n \log{} n)$ computation time, where $n$ is the number of bodies. Canonical use cases include gravity, and charged particle simulations.

This is much faster than a naive N-body approach $O(n^2)$, at high body counts. It is slower than Fast Multiple Methods (FMM), which group target bodies in addition to source ones. $O(n)$ Note that FMM methods are more difficult to implement, but if you have a good implementation, use that instead of this library. There are currently, as of January 2025, no published implementations in Rust.

It uses the [Rayon library](https://docs.rs/rayon/latest/rayon/) to parallelize of the loop summing forces from each source on a given target. You may wish to also parallelize the loop over targets in your application code.

It uses the [lin_alg](https://crates.io/crates/lin_alg) library for the `Vec3` vector type. You will need to import this in your application code, and covert to it when implementing `BodyModel`.

Currently hard-coded for `f64`. Post an issue on GitHub if you'd like `f32` support.

This library is generic over body type and acceleration function. Here's an example of adapting your type for use here:

```rust
use lin_alg::f64::Vec3;

impl BodyModel for Body {
    fn posit(&self) -> Vec3 {
        // Convert to the `lin_alg::f64::Vec3` type here, if not using it directly in your application:
        // Vec3::new(self.posit.x, self.posit.y, self.posit.z)
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
        // Create the tree once per time step.
        let tree = Tree::new(&state.bodies, &bb, &state.config.bh_config);

        // Iterate, in parallel, over target bodies. The loop over source bodies is handled
        // by the acceleration function.
        // bodies.par_iter_mut().enumerate().for_each(|(id, body_target)| { // ...
        for (id, target) in bodies.iter_mut().enumerate() {
            integrate(&config, &tree, target, id);
        }
    }
}

fn integrate(bh_config: &BhConfig, tree: &Tree, target: Body, id_target: usize) {
    // ...

    // This acceleration function can be whatever you'd like. This example shows Newtonian
    // Gravity with MOND.
    let acc_fn = |acc_dir, mass, dist| {
        acc_newton_with_mond(acc_dir, mass, dist, Some(mond_fn), softening_factor_sq)
    };

    let accel = barnes_hut::acc_bh(
        target.posit,
        // `id_target` is used to prevent self-interaction.
        // It is not strictly required, as we also perform a distance check.
        id_target,
        tree,
        bh_config,
        &acc_fn,
    );
}
```