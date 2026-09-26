#!/usr/bin/env bash
set -euo pipefail

python src/s1_comonomers.py
python src/s2_descriptors.py
python src/s3_clean.py
python src/s4_model.py
