use std::cmp::Ordering;
use std::ops::Index;

#[derive(Clone, Debug)]
pub struct Quantity<T> {
    pub value: T,
    pub unit: String,
}

impl<T> Quantity<T> {
    pub fn new(value: T, unit: impl Into<String>) -> Self {
        Self {
            value,
            unit: unit.into(),
        }
    }
}

impl<T> Quantity<Vec<T>> {
    pub fn len(&self) -> usize {
        self.value.len()
    }
}

impl<T: Ord> Quantity<Vec<T>> {
    pub fn sort(&mut self) {
        self.value.sort()
    }
}

impl<T> Quantity<Vec<T>> {
    pub fn sort_by<F>(&mut self, compare: F)
    where
        F: FnMut(&T, &T) -> Ordering,
    {
        self.value.sort_by(compare)
    }
}

impl<T> Index<usize> for Quantity<Vec<T>> {
    type Output = T;

    fn index(&self, index: usize) -> &Self::Output {
        &self.value[index]
    }
}

impl<T: PartialEq> PartialEq for Quantity<T> {
    fn eq(&self, other: &Self) -> bool {
        self.value.eq(&other.value) && self.unit.eq(&other.unit)
    }
}
