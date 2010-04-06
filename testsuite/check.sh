#! /bin/sh
# check script for Lzip - A data compressor based on the LZMA algorithm
# Copyright (C) 2008, 2009, 2010 Antonio Diaz Diaz.
#
# This script is free software: you have unlimited permission
# to copy, distribute and modify it.

LC_ALL=C
export LC_ALL
objdir=`pwd`
testdir=`cd "$1" ; pwd`
LZIP="${objdir}"/lzip
LZIPRECOVER="${objdir}"/lziprecover
framework_failure() { echo "failure in testing framework" ; exit 1 ; }

if [ ! -x "${LZIP}" ] ; then
	echo "${LZIP}: cannot execute"
	exit 1
fi

if [ -d tmp ] ; then rm -rf tmp ; fi
mkdir tmp
printf "testing lzip..."
cd "${objdir}"/tmp

cat "${testdir}"/test1 > in || framework_failure
fail=0

"${LZIP}" -cd "${testdir}"/test1.lz > copy || fail=1
cmp in copy || fail=1

for i in s4Ki 0 1 2 3 4 5 6 7 8 9 ; do
	"${LZIP}" -k -$i in || fail=1
	mv -f in.lz copy.lz || fail=1
	printf "garbage" >> copy.lz || fail=1
	"${LZIP}" -df copy.lz || fail=1
	cmp in copy || fail=1
	printf .
done

for i in s4Ki 0 1 2 3 4 5 6 7 8 9 ; do
	"${LZIP}" -c -$i in > out || fail=1
	printf "g" >> out || fail=1
	"${LZIP}" -cd out > copy || fail=1
	cmp in copy || fail=1
	printf .
done

for i in s4Ki 0 1 2 3 4 5 6 7 8 9 ; do
	"${LZIP}" -$i < in > out || fail=1
	"${LZIP}" -d < out > copy || fail=1
	cmp in copy || fail=1
	printf .
done

for i in s4Ki 0 1 2 3 4 5 6 7 8 9 ; do
	"${LZIP}" -f -$i -o out < in || fail=1
	"${LZIP}" -df -o copy < out.lz || fail=1
	cmp in copy || fail=1
	printf .
done

"${LZIP}" -ce in in in > out || fail=1
printf "garbage" >> out || fail=1
"${LZIPRECOVER}" out || fail=1
for i in 1 2 3 ; do
	"${LZIP}" -cd rec0000${i}out > copy || fail=1
	cmp in copy || fail=1
	printf .
done

echo
if [ ${fail} = 0 ] ; then
	echo "tests completed successfully."
	cd "${objdir}" && rm -r tmp
else
	echo "tests failed."
fi
exit ${fail}
