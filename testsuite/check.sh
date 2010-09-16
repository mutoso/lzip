#! /bin/sh
# check script for Lzip - Data compressor based on the LZMA algorithm
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
if [ ! -x "${LZIPRECOVER}" ] ; then
	echo "${LZIPRECOVER}: cannot execute"
	exit 1
fi

if [ -d tmp ] ; then rm -rf tmp ; fi
mkdir tmp
printf "testing lzip-%s..." "$2"
cd "${objdir}"/tmp

cat "${testdir}"/test.txt > in || framework_failure
fail=0

"${LZIP}" -t "${testdir}"/test_v0.lz || fail=1
printf .
"${LZIP}" -cd "${testdir}"/test_v0.lz > copy || fail=1
cmp in copy || fail=1
printf .

"${LZIP}" -t "${testdir}"/test_v1.lz || fail=1
printf .
"${LZIP}" -cd "${testdir}"/test_v1.lz > copy || fail=1
cmp in copy || fail=1
printf .

"${LZIP}" -t "${testdir}"/test_sync.lz || fail=1
printf .
"${LZIP}" -cd "${testdir}"/test_sync.lz > copy || fail=1
cmp in copy || fail=1
printf .

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

# Description of test files for lziprecover:
# test_bad1.lz: byte at offset 67 changed from 0xCC to 0x33
# test_bad2.lz: [  34-  66) --> copy of bytes [  68- 100)
# test_bad3.lz: [ 512-1536) --> zeroed;       [2560-3584) --> zeroed
# test_bad4.lz: [3072-4096) --> random data;  [4608-5632) --> zeroed
# test_bad5.lz: [1024-2048) --> random data;  [5120-6144) --> random data

printf "\ntesting lziprecover-%s..." "$2"

"${LZIP}" -c in in in > out || fail=1
printf "garbage" >> out || fail=1
"${LZIPRECOVER}" -s out -o out.lz || fail=1
for i in 1 2 3 ; do
	"${LZIP}" -cd rec0000${i}out.lz > copy || fail=1
	cmp in copy || fail=1
	printf .
done

"${LZIP}" -0kf -$i in || fail=1
"${LZIPRECOVER}" -R in.lz > /dev/null || fail=1
printf .
"${LZIPRECOVER}" -R "${testdir}"/test_v1.lz > /dev/null || fail=1
printf .

"${LZIPRECOVER}" -R -o copy.lz "${testdir}"/test_bad1.lz > /dev/null || fail=1
"${LZIP}" -df copy.lz || fail=1
cmp in copy || fail=1
printf .

"${LZIPRECOVER}" -m -o copy.lz "${testdir}"/test_bad1.lz "${testdir}"/test_bad2.lz > /dev/null || fail=1
"${LZIP}" -df copy.lz || fail=1
cmp in copy || fail=1
printf .
"${LZIPRECOVER}" -m -o copy.lz "${testdir}"/test_bad2.lz "${testdir}"/test_bad1.lz > /dev/null || fail=1
"${LZIP}" -df copy.lz || fail=1
cmp in copy || fail=1
printf .

for i in 1 2 ; do
	for j in 3 4 5 ; do
		"${LZIPRECOVER}" -m -o copy.lz "${testdir}"/test_bad${i}.lz "${testdir}"/test_bad${j}.lz > /dev/null || fail=1
		"${LZIP}" -df copy.lz || fail=1
		cmp in copy || fail=1
		printf .
		"${LZIPRECOVER}" -m -o copy.lz "${testdir}"/test_bad${j}.lz "${testdir}"/test_bad${i}.lz > /dev/null || fail=1
		"${LZIP}" -df copy.lz || fail=1
		cmp in copy || fail=1
		printf .
	done
done

"${LZIPRECOVER}" -m -o copy.lz "${testdir}"/test_bad3.lz "${testdir}"/test_bad4.lz "${testdir}"/test_bad5.lz > /dev/null || fail=1
"${LZIP}" -df copy.lz || fail=1
cmp in copy || fail=1
printf .
"${LZIPRECOVER}" -m -o copy.lz "${testdir}"/test_bad4.lz "${testdir}"/test_bad5.lz "${testdir}"/test_bad3.lz > /dev/null || fail=1
"${LZIP}" -df copy.lz || fail=1
cmp in copy || fail=1
printf .
"${LZIPRECOVER}" -m -o copy.lz "${testdir}"/test_bad5.lz "${testdir}"/test_bad3.lz "${testdir}"/test_bad4.lz > /dev/null || fail=1
"${LZIP}" -df copy.lz || fail=1
cmp in copy || fail=1
printf .

echo
if [ ${fail} = 0 ] ; then
	echo "tests completed successfully."
	cd "${objdir}" && rm -r tmp
else
	echo "tests failed."
fi
exit ${fail}
