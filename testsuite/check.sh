#! /bin/sh
# check script for Lzip - LZMA lossless data compressor
# Copyright (C) 2008-2014 Antonio Diaz Diaz.
#
# This script is free software: you have unlimited permission
# to copy, distribute and modify it.

LC_ALL=C
export LC_ALL
objdir=`pwd`
testdir=`cd "$1" ; pwd`
LZIP="${objdir}"/lzip
framework_failure() { echo "failure in testing framework" ; exit 1 ; }

if [ ! -f "${LZIP}" ] || [ ! -x "${LZIP}" ] ; then
	echo "${LZIP}: cannot execute"
	exit 1
fi

if [ -d tmp ] ; then rm -rf tmp ; fi
mkdir tmp
cd "${objdir}"/tmp

cat "${testdir}"/test.txt > in || framework_failure
in_lz="${testdir}"/test.txt.lz
fail=0

printf "testing lzip-%s..." "$2"

"${LZIP}" -cqm4 in > /dev/null
if [ $? = 1 ] ; then printf . ; else printf - ; fail=1 ; fi
"${LZIP}" -cqm274 in > /dev/null
if [ $? = 1 ] ; then printf . ; else printf - ; fail=1 ; fi
"${LZIP}" -cqs-1 in > /dev/null
if [ $? = 1 ] ; then printf . ; else printf - ; fail=1 ; fi
"${LZIP}" -cqs0 in > /dev/null
if [ $? = 1 ] ; then printf . ; else printf - ; fail=1 ; fi
"${LZIP}" -cqs4095 in > /dev/null
if [ $? = 1 ] ; then printf . ; else printf - ; fail=1 ; fi
"${LZIP}" -cqs513MiB in > /dev/null
if [ $? = 1 ] ; then printf . ; else printf - ; fail=1 ; fi
printf "  in: Bad magic number (file not in lzip format).\n" > msg
"${LZIP}" -t in 2> out
if [ $? = 2 ] && cmp out msg ; then printf . ; else printf - ; fail=1 ; fi
printf "  (stdin): Bad magic number (file not in lzip format).\n" > msg
"${LZIP}" -t < in 2> out
if [ $? = 2 ] && cmp out msg ; then printf . ; else printf - ; fail=1 ; fi
rm -f out msg
"${LZIP}" -cdq in
if [ $? = 2 ] ; then printf . ; else printf - ; fail=1 ; fi
"${LZIP}" -cdq < in
if [ $? = 2 ] ; then printf . ; else printf - ; fail=1 ; fi
dd if="${in_lz}" bs=1 count=6 2> /dev/null | "${LZIP}" -tq
if [ $? = 2 ] ; then printf . ; else printf - ; fail=1 ; fi
dd if="${in_lz}" bs=1 count=20 2> /dev/null | "${LZIP}" -tq
if [ $? = 2 ] ; then printf . ; else printf - ; fail=1 ; fi

"${LZIP}" -t "${in_lz}" || fail=1
"${LZIP}" -cd "${in_lz}" > copy || fail=1
cmp in copy || fail=1
printf .

cat "${in_lz}" > copy.lz || framework_failure
printf "to be overwritten" > copy || framework_failure
"${LZIP}" -df copy.lz || fail=1
cmp in copy || fail=1
printf .

printf "to be overwritten" > copy || framework_failure
"${LZIP}" -df -o copy < "${in_lz}" || fail=1
cmp in copy || fail=1
printf .

"${LZIP}" < in > anyothername || fail=1
"${LZIP}" -d anyothername || fail=1
cmp in anyothername.out || fail=1
printf .

cat in in > in2 || framework_failure
"${LZIP}" -o copy2 < in2 || fail=1
"${LZIP}" -t copy2.lz || fail=1
printf .
"${LZIP}" -cd copy2.lz > copy2 || fail=1
cmp in2 copy2 || fail=1
printf .

printf "garbage" >> copy2.lz || framework_failure
printf "to be overwritten" > copy2 || framework_failure
"${LZIP}" -df copy2.lz || fail=1
cmp in2 copy2 || fail=1
printf .

"${LZIP}" -cfq "${in_lz}" > out
if [ $? = 1 ] ; then printf . ; else printf - ; fail=1 ; fi
"${LZIP}" -cF "${in_lz}" > out || fail=1
"${LZIP}" -cd out | "${LZIP}" -d > copy || fail=1
cmp in copy || fail=1
printf .

for i in s4Ki 0 1 2 3 4 5 6 7 8 9 ; do
	"${LZIP}" -k -$i in || fail=1
	mv -f in.lz copy.lz || fail=1
	printf "garbage" >> copy.lz || fail=1
	"${LZIP}" -df copy.lz || fail=1
	cmp in copy || fail=1
done
printf .

for i in s4Ki 0 1 2 3 4 5 6 7 8 9 ; do
	"${LZIP}" -c -$i in > out || fail=1
	printf "g" >> out || fail=1
	"${LZIP}" -cd out > copy || fail=1
	cmp in copy || fail=1
done
printf .

for i in s4Ki 0 1 2 3 4 5 6 7 8 9 ; do
	"${LZIP}" -$i < in > out || fail=1
	"${LZIP}" -d < out > copy || fail=1
	cmp in copy || fail=1
done
printf .

for i in s4Ki 0 1 2 3 4 5 6 7 8 9 ; do
	"${LZIP}" -f -$i -o out < in || fail=1
	"${LZIP}" -df -o copy < out.lz || fail=1
	cmp in copy || fail=1
done
printf .

echo
if [ ${fail} = 0 ] ; then
	echo "tests completed successfully."
	cd "${objdir}" && rm -r tmp
else
	echo "tests failed."
fi
exit ${fail}
