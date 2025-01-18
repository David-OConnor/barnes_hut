//! This algorithm uses Tree Code to group source bodies, as an approximation. It leads to O(N(log N))
//! computation time, where `N` is the number of bodies. Canonical use cases include gravity, and charged
//! particle simulations.
//!
//! See the [readme](https://github.com/David-OConnor/barnes_hut/blob/main/README.md) for details.

use std::{fmt, fmt::Formatter};

#[cfg(feature = "encode")]
use bincode::{Decode, Encode};
use lin_alg::f64::Vec3;
use rayon::prelude::*;

#[derive(Debug)]
#[cfg_attr(feature = "bincode", derive(Encode, Decode))]
pub struct BhConfig {
    /// This determines how aggressively we group. It's no lower than 0. 0 means no grouping.
    /// (Best accuracy; poorest performance; effectively a naive N-body). Higher values
    /// decrease accuracy, and are more performant.
    pub θ: f64,
    pub max_bodies_per_node: usize,
    /// This is a limit on tree division, preventing getting stuck in a loop, e.g. for particles with close.
    /// (or identical) positions
    pub max_tree_depth: usize,
}

impl Default for BhConfig {
    fn default() -> Self {
        Self {
            θ: 0.5,
            max_bodies_per_node: 1,
            max_tree_depth: 15, // todo put back
                                // todo: You have having trouble with the recursion. I think your depth
                                // todo cal logic is causing you to miss sections.
                                // max_tree_depth: 30,
        }
    }
}

/// We use this to allow for arbitrary body (or particle etc) types in application code to
/// use this library.
// If we split this module into a library.
pub trait BodyModel {
    fn posit(&self) -> Vec3;
    fn mass(&self) -> f64;
}

#[derive(Clone, Debug)]
/// A cubical bounding box. length=width=depth.
pub struct Cube {
    pub center: Vec3,
    pub width: f64,
}

impl Cube {
    /// Construct minimum limits that encompass all bodies. Run these each time the bodies change,
    /// or perhaps use a pad and do it at a coarser interval.
    ///
    /// The pad allows us to keep the same cube for multiple timesteps, but taking into acacount
    /// potential movement of bodies outside the cube between these updates.
    ///
    /// The z offset is intended for the case where the Z coordinate for all particles is 0.
    /// This prevents the divisions straddling the points, doubling the number of nodes.
    pub fn from_bodies<T: BodyModel>(bodies: &[T], pad: f64, z_offset: bool) -> Option<Self> {
        if bodies.is_empty() {
            return None;
        }

        let mut x_min = f64::MAX;
        let mut x_max = f64::MIN;
        let mut y_min = f64::MAX;
        let mut y_max = f64::MIN;
        let mut z_min = f64::MAX;
        let mut z_max = f64::MIN;

        for body in bodies {
            let p = &body.posit();
            x_min = x_min.min(p.x);
            x_max = x_max.max(p.x);
            y_min = y_min.min(p.y);
            y_max = y_max.max(p.y);
            z_min = z_min.min(p.z);
            z_max = z_max.max(p.z);
        }

        x_min -= pad;
        x_max += pad;
        y_min -= pad;
        y_max += pad;
        z_min -= pad;
        z_max += pad;

        if z_offset {
            z_max += 1e-5;
        }

        let x_size = x_max - x_min;
        let y_size = y_max - y_min;
        let z_size = z_max - z_min;

        // Coerce to a cube.
        let mut width = x_size.max(y_size).max(z_size);

        let center = Vec3::new(
            (x_max + x_min) / 2.,
            (y_max + y_min) / 2.,
            (z_max + z_min) / 2.,
        );

        Some(Self::new(center, width))
    }

    pub fn new(center: Vec3, width: f64) -> Self {
        Self { center, width }
    }

    /// Divide this into equal-area octants.
    pub(crate) fn divide_into_octants(&self) -> [Self; 8] {
        let width = self.width / 2.;
        let wd2 = self.width / 4.; // short for brevity below.

        // Every combination of + and - for the center offset.
        // The order matters, due to the binary index logic used when partitioning bodies into octants.
        [
            Self::new(self.center + Vec3::new(-wd2, -wd2, -wd2), width),
            Self::new(self.center + Vec3::new(wd2, -wd2, -wd2), width),
            Self::new(self.center + Vec3::new(-wd2, wd2, -wd2), width),
            Self::new(self.center + Vec3::new(wd2, wd2, -wd2), width),
            Self::new(self.center + Vec3::new(-wd2, -wd2, wd2), width),
            Self::new(self.center + Vec3::new(wd2, -wd2, wd2), width),
            Self::new(self.center + Vec3::new(-wd2, wd2, wd2), width),
            Self::new(self.center + Vec3::new(wd2, wd2, wd2), width),
        ]
    }
}

#[derive(Debug)]
pub struct Node {
    /// We use `id` while building the tree, then sort by it, replacing with index.
    /// Once complete, `id` == index in `Tree::nodes`.
    pub id: usize,
    pub bounding_box: Cube,
    /// Node indices in the tree. We use this to guide the transversal process while finding
    /// relevant nodes for a given target body.
    pub children: Vec<usize>,
    mass: f64,
    center_of_mass: Vec3,
}

impl fmt::Display for Node {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "Id: {}, Width: {:.3}, Ch: {:?}",
            self.id, self.bounding_box.width, self.children
        )
    }
}

#[derive(Debug)]
/// A recursive tree. Each node can be subdivided  Terminates with `NodeType::NodeTerminal`.
pub struct Tree {
    // pub struct Tree<'a> {
    /// Order matters; we index this by `Node::children`.
    pub nodes: Vec<Node>,
    // pub nodes: &'a mut Vec<Node>,
}

impl Tree {
    /// Constructs a tree. Call this externaly using all bodies, once per time step.
    /// It creates the entire tree, branching until each cell has `MAX_BODIES_PER_NODE` or fewer
    /// bodies, or it reaches a maximum recursion depth.
    ///
    /// We partially transverse it as-required while calculating the force on a given target.
    // pub fn new(bodies: &[Body], bb: &Cube, config: &BhConfig) -> Self {
    pub fn new<T: BodyModel>(bodies: &[T], bb: &Cube, config: &BhConfig) -> Self {
        // Convert &[Body] to &[&Body].
        // let body_refs: Vec<&Body> = bodies.iter().collect();
        let body_refs: Vec<&T> = bodies.iter().collect();

        // todo: Refine this guess A/R.
        // From an unrigorous benchmark, preallocating seems to be slightly faster, but not significantly so?
        let mut nodes = Vec::with_capacity(bodies.len() * 7 / 4);

        let mut current_node_i: usize = 0;

        // Stack to simulate recursion: Each entry contains (bodies, bounding box, parent_id, child_index, depth).
        let mut stack = Vec::new();
        stack.push((body_refs.to_vec(), bb.clone(), None, 0));

        while let Some((bodies_, bb_, parent_id, depth)) = stack.pop() {
            if depth > config.max_tree_depth {
                // eprintln!("Tree generation: Max recusion depth exceeded.");
                break;
            }
            let (center_of_mass, mass) = center_of_mass(&bodies_);

            let node_id = current_node_i;
            nodes.push(Node {
                id: node_id,
                bounding_box: bb_.clone(),
                mass,
                center_of_mass,
                children: Vec::new(),
            });

            current_node_i += 1;

            if let Some(pid) = parent_id {
                // Rust is requesting an explicit type here.
                let n: &mut Node = &mut nodes[pid];
                n.children.push(node_id);
            }

            // If multiple (past our threshold) bodies are in this node, create an internal node and push its ID.
            // Divide into octants and partition bodies. Otherwise, create a leaf node.
            if bodies_.len() > config.max_bodies_per_node {
                let octants = bb_.divide_into_octants();
                let bodies_by_octant = partition(&bodies_, &bb_);

                // Add each octant with bodies to the stack.
                for (i, octant) in octants.into_iter().enumerate() {
                    if !bodies_by_octant[i].is_empty() {
                        stack.push((
                            bodies_by_octant[i].clone(),
                            octant,
                            Some(node_id),
                            depth + 1,
                        ));
                    }
                }
            }
        }

        // Now that nodes are populated, rearrange so index == `id`. We will then index by `children`.
        nodes.sort_by(|l, r| l.id.partial_cmp(&r.id).unwrap());

        Self { nodes }
    }

    /// Get all leaves relevant to a given target. We use this to create a coarser
    /// version of the tree, containing only the nodes we need to calculate acceleration
    /// on a specific target.
    pub fn leaves(&self, posit_target: Vec3, id_target: usize, config: &BhConfig) -> Vec<&Node> {
        let mut result = Vec::new();

        if self.nodes.is_empty() {
            return result;
        }

        let node_i = 0;

        let mut stack = Vec::new();
        stack.push(node_i);

        while let Some(current_node_i) = stack.pop() {
            let node = &self.nodes[current_node_i];

            if node.children.len() <= config.max_bodies_per_node {
                result.push(node);
                continue;
            }

            let dist = (posit_target - node.center_of_mass).magnitude();

            // Avoid self-interaction based on distance or id_target.
            // todo: Use id_target, if able.
            if dist < 1e-10 {
                continue;
            }

            if node.bounding_box.width / dist < config.θ {
                result.push(node);
            } else {
                // The source is near; add children to the stack to go deeper.
                for &child_i in &node.children {
                    stack.push(child_i);
                }
            }
        }

        result
    }
}

/// Compute center of mass as a position, and mass value.
// fn center_of_mass(bodies: &[&Body]) -> (Vec3, f64) {
fn center_of_mass<T: BodyModel>(bodies: &[&T]) -> (Vec3, f64) {
    let mut mass = 0.;
    let mut center_of_mass = Vec3::new_zero();

    for body in bodies {
        mass += body.mass();
        // Weight the center by mass.
        center_of_mass += body.posit() * body.mass();
    }
    if mass.abs() > f64::EPSILON {
        // Remove the weighting's effect on magnitude of the center.
        center_of_mass /= mass;
    }

    (center_of_mass, mass)
}

/// Partition bodies into each of the 8 octants.
// fn partition<'a>(bodies: &[&'a Body], bb: &Cube) -> [Vec<&'a Body>; 8] {
fn partition<'a, T: BodyModel>(bodies: &[&'a T], bb: &Cube) -> [Vec<&'a T>; 8] {
    let mut result: [Vec<&T>; 8] = Default::default();

    for body in bodies {
        let mut index = 0;
        if body.posit().x > bb.center.x {
            index |= 0b001;
        }
        if body.posit().y > bb.center.y {
            index |= 0b010;
        }
        if body.posit().z > bb.center.z {
            index |= 0b100;
        }

        result[index].push(body);
    }

    result
}

type AccFn<'a> = dyn Fn(Vec3, f64, f64) -> Vec3 + Send + Sync + 'a;

/// Calculate Newtonian acceleration using the Barnes Hut algorithm.
pub fn acc_bh(
    posit_target: Vec3,
    id_target: usize,
    tree: &Tree,
    config: &BhConfig,
    acc_fn: &AccFn,
) -> Vec3 {
    // todo: Put back the part checking for self interaction.
    tree.leaves(posit_target, id_target, config)
        .par_iter()
        .filter_map(|leaf| {
            let acc_diff = leaf.center_of_mass - posit_target;
            let dist = acc_diff.magnitude();

            // todo: Not sure why we get 0 here when we check it in `leaves()`.
            // todo: QC this.
            if dist < 1e-8 {
                return None;
            }

            let acc_dir = acc_diff / dist; // Unit vec

            Some(acc_fn(acc_dir, leaf.mass, dist))
        })
        .reduce(Vec3::new_zero, |acc, elem| acc + elem)
}
