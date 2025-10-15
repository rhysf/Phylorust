#[derive(Clone)]
pub struct Logger;

use colored::*;

impl Logger {
    pub fn output(&self, message : &str) {
        println!("{}", message);
    }

    pub fn information(&self, message : &str) {
        eprintln!("{}", message);
    }

    pub fn warning(&self, message : &str) {
        eprintln!("{}", message.yellow());
    }

    pub fn error(&self, message : &str) {
        eprintln!("{}", message.red());
    }
}