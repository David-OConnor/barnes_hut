# Runs the Barnes Hut algorithm for fast n-body simulations.

[![Crate](https://img.shields.io/crates/v/barnes_hut.svg)](https://crates.io/crates/barnes_hut)
[![Docs](https://docs.rs/lin_alg/badge.svg)](https://docs.rs/barnes_hut)

This algorithm uses Tree Code to group source bodies, as an approximation. It leads to O(N(log N)) computation time, where `N` is the number of bodies. Canonical use cases include gravity, and charged particle simulations.

This is much faster than a naive N-body approach, at high body counts. It is slower than Fast Multiple Methods (FMM),
which group target bodies in addition to source ones. Note that FMM methods are more difficult to implement, but if you have a 
good implementation, use that instead of this library. There are currently, as of January 2025, no published implementations
in Rust.

Uses the Rayon library for parallelization. 

Currently hard-coded for `f64`. Post an issue on Github if you'd like `f32` support.

This library is generic over `Body` and acceleration function . Here's an example of adapting your type for use here:

```rust
impl BodyModel for Body {
    fn posit(&self) -> Vec3 {
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
    
    let tree = Tree::new(&state.bodies, &bb, &state.config.bh_config);
    
    // ...
    integrate(&config, &tree);
}

fn integrate(bh_config: &BhConfig, tree: &Tree) {
    // ...

    // This acceleration function can be whatever you'd like. This example shows Newtonian
    // Gravity with MOND.
    let acc_fn = |acc_dir, mass, dist| {
        acc_newton_with_mond(acc_dir, mass, dist, Some(mond_fn), softening_factor_sq)
    };

    let accel = barnes_hut::acc_bh(
        posit_target,
        id_target, // Currently unused; set to 0 if you'd like, or the enumeration of bodies.
        tree,
        bh_config,
        Box::new(acc_fn),
    );
}
```