//! This algorithm uses Tree Code to group source bodies, as an approximation. It leads to O(N(log N))
//! computation time, where `N` is the number of bodies. Canonical use cases include gravity, and charged
//! particle simulations.
//!
//! See the [readme](https://github.com/David-OConnor/barnes_hut/blob/main/README.md) for details,
//! including an example.

#![allow(non_ascii_idents)]
#![allow(mixed_script_confusables)]

// todo: Ideally make generic over f32 and f64, but we don't have a good way to do that with `Vec3`.

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
/// use this library. Substitute `charge` for `mass` as required.
pub trait BodyModel {
    fn posit(&self) -> Vec3;
    fn mass(&self) -> f64;
}

#[derive(Clone, Debug)]
#[cfg_attr(feature = "bincode", derive(Encode, Decode))]
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
        let width = x_size.max(y_size).max(z_size);

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
    /// Mass, center-of-mass, and body_ids include those from all sub-nodes.
    pub id: usize,
    pub bounding_box: Cube,
    /// Node indices in the tree. We use this to guide the transversal process while finding
    /// relevant nodes for a given target body.
    pub children: Vec<usize>,
    pub mass: f64,
    pub center_of_mass: Vec3,
    pub body_ids: Vec<usize>,
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
    /// Order matters; we index this by `Node::children`.
    // Note: It doesn't appear that passing in a persistent, pre-allocated nodes Vec from the applicatoni
    // has a significant impact on tree construction time.
    pub nodes: Vec<Node>,
}

impl Tree {
    /// Constructs a tree. Call this externaly using all bodies, once per time step.
    /// It creates the entire tree, branching until each cell has `MAX_BODIES_PER_NODE` or fewer
    /// bodies, or it reaches a maximum recursion depth.
    ///
    /// We partially transverse it as-required while calculating the force on a given target.
    pub fn new<T: BodyModel>(bodies: &[T], bb: &Cube, config: &BhConfig) -> Self {
        // Convert &[T] to &[&T].
        let body_refs: Vec<&T> = bodies.iter().collect();

        // todo: Refine this guess A/R.
        // From an unrigorous benchmark, preallocating seems to be slightly faster, but not significantly so?
        let mut nodes = Vec::with_capacity(bodies.len() * 7 / 4);

        let mut current_node_i: usize = 0;

        // Stack to simulate recursion: Each entry contains (bodies, bounding box, parent_id, child_index, depth).
        let mut stack = Vec::new();

        // body ids matches indexes with bodies.
        let body_ids_init: Vec<usize> = body_refs.iter().enumerate().map(|(id, _)| id).collect();

        stack.push((body_refs.to_vec(), body_ids_init, bb.clone(), None, 0));

        while let Some((bodies_, body_ids, bb_, parent_id, depth)) = stack.pop() {
            if depth > config.max_tree_depth {
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
                body_ids: body_ids.clone(), // todo: The clone...
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
                let bodies_by_octant = partition(&bodies_, &body_ids, &bb_);

                // Add each octant with bodies to the stack.
                for (i, octant) in octants.into_iter().enumerate() {
                    if !bodies_by_octant[i].is_empty() {
                        let mut bto = Vec::with_capacity(bodies_by_octant[i].len());
                        let mut ids_this_octant = Vec::with_capacity(bodies_by_octant[i].len());

                        // todo: The clone etc?
                        for (body, id) in &bodies_by_octant[i] {
                            bto.push(*body);
                            ids_this_octant.push(*id);
                        }

                        stack.push((bto, ids_this_octant, octant, Some(node_id), depth + 1));
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
    pub fn leaves(&self, posit_target: Vec3, config: &BhConfig) -> Vec<&Node> {
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
fn partition<'a, T: BodyModel>(
    bodies: &[&'a T],
    body_ids: &[usize],
    bb: &Cube,
) -> [Vec<(&'a T, usize)>; 8] {
    let mut result: [Vec<(&'a T, usize)>; 8] = Default::default();

    for (i, body) in bodies.iter().enumerate() {
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

        result[index].push((body, body_ids[i]));
    }

    result
}

/// Calculate force using the Barnes Hut algorithm. The force function passed
/// as a parameter has signature `(acc_dir: Vec3 (unit), mass_src: f64, distance: f64) -> Vec3`
/// `id_target` is the index in the body array used to make the tree; it prevents self-interaction.
/// Note that `mass` can be interchanged with `charge`, or similar.
///
/// When handling target mass or charge, reflect that in your `force_fn`; not here.
pub fn run_bh<F>(
    posit_target: Vec3,
    id_target: usize,
    tree: &Tree,
    config: &BhConfig,
    force_fn: &F,
) -> Vec3
where
    F: Fn(Vec3, f64, f64) -> Vec3 + Send + Sync,
{
    tree.leaves(posit_target, config)
        .par_iter()
        .filter_map(|leaf| {
            if leaf.body_ids.contains(&id_target) {
                // Prevent self-interaction.
                return None;
            }

            let acc_diff = leaf.center_of_mass - posit_target;
            let dist = acc_diff.magnitude();

            let acc_dir = acc_diff / dist; // Unit vec

            Some(force_fn(acc_dir, leaf.mass, dist))
        })
        .reduce(Vec3::new_zero, |acc, elem| acc + elem)
}
