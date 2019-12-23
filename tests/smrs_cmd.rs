use std::fs;
use std::process::Command;

use assert_cmd::prelude::*;
use predicates::str::contains;
use tempfile::TempDir;

#[test]
fn search() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("smrs")?;

    cmd.arg("search")
        .arg("tests/test-data/demo/SRR2060939_1.sig")
        .arg("tests/test-data/v5.sbt.json")
        .assert()
        .success()
        .stdout(contains("SRR2060939_1.fastq.gz"))
        .stdout(contains("SRR2060939_2.fastq.gz"))
        .stdout(contains("SRR2255622_1.fastq.gz"));

    Ok(())
}

#[test]
#[ignore]
fn search_only_leaves() -> Result<(), Box<dyn std::error::Error>> {
    let mut cmd = Command::cargo_bin("smrs")?;

    cmd.arg("search")
        .arg("tests/test-data/demo/SRR2060939_1.sig")
        .arg("tests/test-data/leaves.sbt.json")
        .assert()
        .success()
        .stdout(contains("SRR2060939_1.fastq.gz"))
        .stdout(contains("SRR2060939_2.fastq.gz"))
        .stdout(contains("SRR2255622_1.fastq.gz"));

    Ok(())
}

#[test]
#[ignore]
#[cfg(unix)]
fn compute_index_and_search() -> Result<(), Box<dyn std::error::Error>> {
    let tmp_dir = TempDir::new()?;
    fs::copy("tests/test-data/short.fa", tmp_dir.path().join("short.fa"))?;
    fs::copy(
        "tests/test-data/short2.fa",
        tmp_dir.path().join("short2.fa"),
    )?;

    assert!(tmp_dir.path().join("short.fa").exists());
    assert!(tmp_dir.path().join("short2.fa").exists());

    let mut cmd = Command::new("sourmash");
    cmd.arg("compute")
        .args(&["short.fa", "short2.fa"])
        .current_dir(&tmp_dir)
        .assert()
        .success();

    assert!(tmp_dir.path().join("short.fa.sig").exists());
    assert!(tmp_dir.path().join("short2.fa.sig").exists());

    let mut cmd = Command::new("sourmash");
    //let mut cmd = Command::cargo_bin("smrs")?;
    cmd.arg("index")
        .args(&["-k", "31"])
        //.args(&["-o", "zzz.sbt.json"])
        .arg("zzz.sbt.json")
        .args(&["short.fa.sig", "short2.fa.sig"])
        .current_dir(&tmp_dir)
        .assert()
        .success();

    assert!(tmp_dir.path().join("zzz.sbt.json").exists());

    let cmds = vec![Command::new("sourmash"), Command::cargo_bin("smrs")?];

    for mut cmd in cmds {
        cmd.arg("search")
            .args(&["-k", "31"])
            .arg("short.fa.sig")
            .arg("zzz.sbt.json")
            .current_dir(&tmp_dir)
            .assert()
            .success()
            .stdout(contains("short.fa"))
            .stdout(contains("short2.fa"));
    }

    Ok(())
}

#[test]
#[cfg(unix)]
fn index_and_search() -> Result<(), Box<dyn std::error::Error>> {
    let tmp_dir = TempDir::new()?;
    fs::copy(
        "tests/test-data/demo/SRR2060939_1.sig",
        tmp_dir.path().join("1.sig"),
    )?;
    fs::copy(
        "tests/test-data/demo/SRR2060939_2.sig",
        tmp_dir.path().join("2.sig"),
    )?;

    assert!(tmp_dir.path().join("1.sig").exists());
    assert!(tmp_dir.path().join("2.sig").exists());

    let mut cmd = Command::cargo_bin("smrs")?;
    cmd.arg("index")
        .args(&["-k", "31"])
        .args(&["-o", "zzz.sbt.json"])
        .args(&["1.sig", "2.sig"])
        .current_dir(&tmp_dir)
        .assert()
        .success();

    assert!(tmp_dir.path().join("zzz.sbt.json").exists());

    let cmds = vec![Command::new("sourmash"), Command::cargo_bin("smrs")?];

    for mut cmd in cmds {
        cmd.arg("search")
            .args(&["-k", "31"])
            .arg("1.sig")
            .arg("zzz.sbt.json")
            .current_dir(&tmp_dir)
            .assert()
            .success()
            .stdout(contains("2 matches:"))
            .stdout(contains("SRR2060939_1.fastq.gz"))
            .stdout(contains("SRR2060939_2.fastq.gz"));
    }

    Ok(())
}
