#!/bin/sh

# Script to stage the document files (generated in place) but not the sources
find . \
  -mindepth 1 \
  -type f \( -name 'CMakeLists.txt' -prune -o -name 'stage_doc_build.sh' \) \
  -o -type d \( -path './doxygen' -prune -o -path './images' -prune -o -path './sphinx' -prune \) \
  -o -print \
| xargs -i git add {}
