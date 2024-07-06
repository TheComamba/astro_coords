#!/bin/bash

cargo install cargo-llvm-cov
sudo apt-get update
sudo apt-get install lcov -y

# Doctests require nightly and are not included in --all-targets
# Compare: https://github.com/xd009642/tarpaulin/issues/850
rustup update nightly
time cargo +nightly llvm-cov --all-features --workspace --lcov --doctests --output-path lcov.info
genhtml lcov.info --output-directory code_coverage
