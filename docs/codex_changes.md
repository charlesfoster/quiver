# Codex Changes

## Docker Mac LoFreq Stability

Added a new `docker_mac` profile for Docker Desktop on Apple Silicon. It inherits the standard Docker profile, then sets `params.lofreq_pp_threads = 1` and constrains `LOFREQ_CALL` to one CPU and one concurrent fork. This makes the profile use serial `lofreq call` instead of LoFreq's `call-parallel` wrapper, avoiding the repeated Docker/Rosetta SIGKILL failures seen when multiple amd64 LoFreq child processes run at once.

Added `params.lofreq_pp_threads` as a configurable pipeline parameter. The default remains `8`, preserving the previous behavior on non-Mac profiles. When the value is greater than one, `LOFREQ_CALL` runs `lofreq call-parallel` with the requested worker count capped by `task.cpus`. When the value is one, it runs plain serial `lofreq call`.

Removed the `lofreq alnqual` preprocessing integration. `LOFREQ_PREPROCESS` now runs only `lofreq indelqual --dindel`, indexes that BAM, and passes the indel-quality BAM downstream. This removes another LoFreq subcommand that was being killed under Docker Desktop on Apple Silicon while keeping the BI/BD indel-quality tags required for `--call-indels`.
