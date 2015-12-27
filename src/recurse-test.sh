#!/bin/sh
set -e
# Check that the output of test_parse matches expectations
# Also that that the output of test_parse can be re-parsed
# and give the same output.

prog="$1"
shift
input="$1"
shift
expect="$1"
shift

ibase="$(basename "$input")"

trap 'rm -f "$ibase".out1 "$ibase".out2 "$ibase".err1 "$ibase".err2' INT TERM QUIT EXIT

"$prog" "$input" > "$ibase".out1 2> "$ibase".err1

"$prog" "$ibase".out1 > "$ibase".out2 2> "$ibase".err2

ret=0

echo "Difference expected -> pass 1"

if ! diff -u "$expect" "$ibase".out1
then
  ret=1
fi

echo "Difference pass 1 -> pass 2"

if ! diff -u "$ibase".out1 "$ibase".out2
then
  ret=1
fi

if [ $ret -ne 0 ]; then
  echo "============ Pass 1 stderr ============"
  tail -n100 "$ibase".err1
  echo "============ Pass 2 stderr ============"
  tail -n100 "$ibase".err2
fi

exit $ret
