#!/usr/bin/env bash
set -euo pipefail
set -x
# -----------------------------------------------------------------------------
# Driver script to compile and run main_const_pulse_multi_piece_eem3
# Requires: g++, source file in cwd, dlib headers
# -----------------------------------------------------------------------------

# log() {
#   echo "[$(date +'%Y-%m-%d %H:%M:%S')] $*"
# }

# Configuration
SRC="main_const_pulse_multi_piece.cpp"
EXE="main_const_pulse_multi_piece.exe"
OUT_DIR="new_runs_log"
ampof_values=(0.035 )
q_values=(12 13 14 15)
target_values=(19 20)
THREADS=32

# # Preconditions
# command -v g++ >/dev/null || { log "Error: g++ not installed."; exit 1; }
# [[ -f "$SRC" ]] || { log "Error: Source file '$SRC' not found."; exit 1; }

# # Compile
# log "Compiling $SRC..."
# g++ "$SRC" -Idlib/dlib -std=c++17 -O2 -o "$EXE"
# log "Compilation successful."

# Run all combinations
for ampof in "${ampof_values[@]}"; do
  for q in "${q_values[@]}"; do
    for max_target in "${target_values[@]}"; do
      outfile="${OUT_DIR}/test_amp_${ampof}_q${q}_n${max_target}_p${THREADS}.txt"
      log "Running: ampof=$ampof, q=$q, max_target=$max_target"
      "./$EXE" "$ampof" "$q" "$max_target" "$THREADS" > "$outfile"
    done
  done
done

log "All runs completed."

echo "Press ENTER to close..."
read -r
