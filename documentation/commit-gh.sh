#!/bin/sh
set -e
# Usage: commit-gh <repo/.git> <files>
#
# Creates a commit containing only the files provided as arguments
#
# Does not disturb the working copy or index

export GIT_DIR="$1"
shift

# Commit to this branch
BRANCH=refs/heads/gh-pages

# Use the main branch description as the gh-pages commit message
MSG=`git describe --tags --always`

# Scratch space
TDIR=`mktemp -d -p $PWD`

# Automatic cleanup of scratch space
trap 'rm -rf $TDIR' INT TERM QUIT EXIT

export GIT_INDEX_FILE="$TDIR/index"

# Add listed files to a new (empty) index
git update-index --add "$@"

# Hash and store a zero length blob
NJ=`git hash-object --stdin -w < /dev/null`

# add the magic (empty) file which tells github not to run its html generator
git update-index --add --cacheinfo 100644,$NJ,.nojekyll

# Write the index into the repo, get tree hash
TREE=`git write-tree`

echo "TREE $TREE"
git cat-file -p $TREE

# Create a commit with our new tree
CMT=`git commit-tree -m "$MSG" $TREE`

echo "COMMIT $CMT"
git cat-file -p $CMT

# Update the branch with the new commit tree hash
git update-ref $BRANCH $CMT

echo "Done"
