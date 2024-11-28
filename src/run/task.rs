use serde::{Serialize, Deserialize};

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Task {
    /// A unique identifier for the task (optional, useful for debugging or logging)
    pub id: usize,
    /// Whether this task should permute the data (true for permutations)
    pub permute: bool,
}
