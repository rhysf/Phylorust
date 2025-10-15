use std::sync::{Arc, Mutex};

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

// Extension trait to forward calls through Arc<Mutex<Logger>>
pub trait LogExt {
    fn information(&self, msg: &str);
    fn warning(&self, msg: &str);
    fn error(&self, msg: &str);
    //fn output(&self, msg: &str);
}

// Support both Arc<Mutex<Logger>> and &Arc<Mutex<Logger>>
impl LogExt for Arc<Mutex<Logger>> {
    fn information(&self, msg: &str) { self.lock().unwrap().information(msg); }
    fn warning(&self, msg: &str)     { self.lock().unwrap().warning(msg); }
    fn error(&self, msg: &str)       { self.lock().unwrap().error(msg); }
    //fn output(&self, msg: &str)      { self.lock().unwrap().output(msg); }
}

impl LogExt for &Arc<Mutex<Logger>> {
    fn information(&self, msg: &str) { self.lock().unwrap().information(msg); }
    fn warning(&self, msg: &str)     { self.lock().unwrap().warning(msg); }
    fn error(&self, msg: &str)       { self.lock().unwrap().error(msg); }
    //fn output(&self, msg: &str)      { self.lock().unwrap().output(msg); }
}