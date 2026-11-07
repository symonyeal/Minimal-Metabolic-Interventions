//! Write-once result files. Stage and sync complete bytes, then publish without replacing a name.
//! A solution file is published before its summary; only the summary accepts the pair.

use std::fs::{self, File};
use std::io::{self, Read, Write};
use std::path::{Path, PathBuf};
use std::time::{SystemTime, UNIX_EPOCH};

use sha2::{Digest, Sha256};

pub struct Pending {
    pub file: File,
    pub dest: PathBuf,
    path: PathBuf,
}

impl Pending {
    pub fn new(path: &str) -> io::Result<Self> {
        let path = Path::new(path);
        let parent = path.parent().filter(|p| !p.as_os_str().is_empty()).unwrap_or(Path::new("."));
        let dest =
            parent.canonicalize()?.join(path.file_name().ok_or_else(|| io::Error::other("output needs a filename"))?);
        if dest.exists() {
            return Err(io::Error::new(
                io::ErrorKind::AlreadyExists,
                format!("{} exists; choose a new result filename", dest.display()),
            ));
        }
        let id = SystemTime::now().duration_since(UNIX_EPOCH).unwrap().as_nanos();
        let path = parent.join(format!(".mcs-{}-{id}.pending", std::process::id()));
        let file = File::options().read(true).write(true).create_new(true).open(&path)?;
        Ok(Self { file, dest, path })
    }

    pub fn digest(&mut self) -> io::Result<String> {
        use std::io::{Seek, SeekFrom};
        self.file.flush()?;
        self.file.sync_all()?;
        self.file.seek(SeekFrom::Start(0))?;
        let mut hash = Sha256::new();
        let mut b = [0; 64 * 1024];
        loop {
            let n = self.file.read(&mut b)?;
            if n == 0 {
                break;
            }
            hash.update(&b[..n]);
        }
        Ok(hash.finalize().iter().map(|b| format!("{b:02x}")).collect())
    }

    pub fn publish(self) -> io::Result<()> {
        self.file.sync_all()?;
        // An atomic create-if-absent: another writer can never replace an accepted result.
        fs::hard_link(&self.path, &self.dest)
    }
}

impl Drop for Pending {
    fn drop(&mut self) {
        let _ = fs::remove_file(&self.path);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn two_staged_writers_cannot_replace_one_result() {
        let root = std::env::var_os("MCS_TEST_ROOT")
            .map(PathBuf::from)
            .unwrap_or_else(|| Path::new(env!("CARGO_MANIFEST_DIR")).join("target/output-tests"));
        let dir = root.join(format!(
            "{}-{}",
            std::process::id(),
            SystemTime::now().duration_since(UNIX_EPOCH).unwrap().as_nanos()
        ));
        fs::create_dir_all(&dir).unwrap();
        let dest = dir.join("result.json");
        let mut a = Pending::new(dest.to_str().unwrap()).unwrap();
        let mut b = Pending::new(dest.to_str().unwrap()).unwrap();
        a.file.write_all(b"accepted").unwrap();
        b.file.write_all(b"other").unwrap();
        a.publish().unwrap();
        assert_eq!(b.publish().unwrap_err().kind(), io::ErrorKind::AlreadyExists);
        assert_eq!(fs::read(dest).unwrap(), b"accepted");
    }
}
