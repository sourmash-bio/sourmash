use std::cell::RefCell;
use std::ffi::CStr;
use std::mem;
use std::os::raw::c_char;
use std::panic;
use std::ptr;
use std::slice;
use std::str;
use std::thread;

use crate::errors::SourmashErrorCode;
use failure::{Error, Fail};

thread_local! {
    pub static LAST_ERROR: RefCell<Option<Error>> = RefCell::new(None);
}

macro_rules! ffi_fn (
    // a function that catches panics and returns a result (err goes to tls)
    (
        $(#[$attr:meta])*
        unsafe fn $name:ident($($aname:ident: $aty:ty),* $(,)*) -> Result<$rv:ty> $body:block
    ) => (
        #[no_mangle]
        $(#[$attr])*
        pub unsafe extern "C" fn $name($($aname: $aty,)*) -> $rv
        {
            $crate::ffi::utils::landingpad(|| $body)
        }
    );

    // a function that catches panics and returns nothing (err goes to tls)
    (
        $(#[$attr:meta])*
        unsafe fn $name:ident($($aname:ident: $aty:ty),* $(,)*) $body:block
    ) => {
        #[no_mangle]
        $(#[$attr])*
        pub unsafe extern "C" fn $name($($aname: $aty,)*)
        {
            // this silences panics and stuff
            $crate::ffi::utils::landingpad(|| { $body; Ok(0 as ::std::os::raw::c_int) });
        }
    }
);

/// An error thrown by `landingpad` in place of panics.
#[derive(Fail, Debug)]
#[fail(display = "sourmash panicked: {}", _0)]
pub struct Panic(String);

/// Returns the last error message.
///
/// If there is no error an empty string is returned.  This allocates new memory
/// that needs to be freed with `sourmash_str_free`.
#[no_mangle]
pub unsafe extern "C" fn sourmash_err_get_last_message() -> SourmashStr {
    use std::fmt::Write;
    LAST_ERROR.with(|e| {
        if let Some(ref err) = *e.borrow() {
            let mut msg = err.to_string();
            for cause in err.iter_causes() {
                write!(&mut msg, "\n  caused by: {}", cause).ok();
            }
            SourmashStr::from_string(msg)
        } else {
            Default::default()
        }
    })
}

/// Returns the panic information as string.
#[no_mangle]
pub unsafe extern "C" fn sourmash_err_get_backtrace() -> SourmashStr {
    LAST_ERROR.with(|e| {
        if let Some(ref error) = *e.borrow() {
            let backtrace = error.backtrace().to_string();
            if !backtrace.is_empty() {
                use std::fmt::Write;
                let mut out = String::new();
                write!(&mut out, "stacktrace: {}", backtrace).ok();
                SourmashStr::from_string(out)
            } else {
                Default::default()
            }
        } else {
            Default::default()
        }
    })
}

/// Clears the last error.
#[no_mangle]
pub unsafe extern "C" fn sourmash_err_clear() {
    LAST_ERROR.with(|e| {
        *e.borrow_mut() = None;
    });
}

/// Initializes the library
#[no_mangle]
pub unsafe extern "C" fn sourmash_init() {
    set_panic_hook();
}

/// Returns the last error code.
///
/// If there is no error, 0 is returned.
#[no_mangle]
pub unsafe extern "C" fn sourmash_err_get_last_code() -> SourmashErrorCode {
    LAST_ERROR.with(|e| {
        if let Some(ref err) = *e.borrow() {
            SourmashErrorCode::from_error(err)
        } else {
            SourmashErrorCode::NoError
        }
    })
}

fn set_last_error(err: Error) {
    LAST_ERROR.with(|e| {
        *e.borrow_mut() = Some(err);
    });
}

pub unsafe fn set_panic_hook() {
    panic::set_hook(Box::new(|info| {
        let thread = thread::current();
        let thread = thread.name().unwrap_or("unnamed");

        let message = match info.payload().downcast_ref::<&str>() {
            Some(s) => *s,
            None => match info.payload().downcast_ref::<String>() {
                Some(s) => &**s,
                None => "Box<Any>",
            },
        };

        let description = match info.location() {
            Some(location) => format!(
                "thread '{}' panicked with '{}' at {}:{}",
                thread,
                message,
                location.file(),
                location.line()
            ),
            None => format!("thread '{}' panicked with '{}'", thread, message),
        };

        set_last_error(Panic(description).into())
    }));
}

pub unsafe fn landingpad<F, T>(f: F) -> T
where
    F: FnOnce() -> Result<T, Error> + panic::UnwindSafe,
{
    match panic::catch_unwind(f) {
        Ok(Ok(result)) => result,
        Ok(Err(err)) => {
            set_last_error(err);
            mem::zeroed()
        }
        Err(_) => mem::zeroed(),
    }
}

/// Represents a string.
#[repr(C)]
pub struct SourmashStr {
    pub data: *mut c_char,
    pub len: usize,
    pub owned: bool,
}

impl Default for SourmashStr {
    fn default() -> SourmashStr {
        SourmashStr {
            data: ptr::null_mut(),
            len: 0,
            owned: false,
        }
    }
}

impl SourmashStr {
    pub fn new(s: &str) -> SourmashStr {
        SourmashStr {
            data: s.as_ptr() as *mut c_char,
            len: s.len(),
            owned: false,
        }
    }

    pub fn from_string(mut s: String) -> SourmashStr {
        s.shrink_to_fit();
        let rv = SourmashStr {
            data: s.as_ptr() as *mut c_char,
            len: s.len(),
            owned: true,
        };
        Box::into_raw(s.into_boxed_str());
        rv
    }

    pub unsafe fn free(&mut self) {
        if self.owned {
            String::from_raw_parts(self.data as *mut _, self.len, self.len);
            self.data = ptr::null_mut();
            self.len = 0;
            self.owned = false;
        }
    }

    pub fn as_str(&self) -> &str {
        unsafe { str::from_utf8_unchecked(slice::from_raw_parts(self.data as *const _, self.len)) }
    }
}

ffi_fn! {
    /// Creates a sourmash str from a c string.
    ///
    /// This sets the string to owned.  In case it's not owned you either have
    /// to make sure you are not freeing the memory or you need to set the
    /// owned flag to false.
    unsafe fn sourmash_str_from_cstr(s: *const c_char) -> Result<SourmashStr> {
        let s = CStr::from_ptr(s).to_str()?;
        Ok(SourmashStr {
            data: s.as_ptr() as *mut _,
            len: s.len(),
            owned: true,
        })
    }
}

/// Frees a sourmash str.
///
/// If the string is marked as not owned then this function does not
/// do anything.
#[no_mangle]
pub unsafe extern "C" fn sourmash_str_free(s: *mut SourmashStr) {
    if !s.is_null() {
        (*s).free()
    }
}
