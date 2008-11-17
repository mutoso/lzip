/*  Unzcrash - A test program written to test robustness to
               decompression of corrupted data.
    Inspired by unzcrash.c from Julian Seward's bzip2.
    Copyright (C) 2008 Antonio Diaz Diaz.

    This program is free software: you have unlimited permission
    to copy, distribute and modify it.

    Usage is:
      unzcrash "lzip -tv" filename.lz

    This program reads the specified file and then repeatedly
    decompresses it, each time with a different bit of the compressed
    data inverted, so as to test all possible one-bit errors. This
    should not cause any invalid memory accesses. If it does, please,
    report it as a bug.

    Compile this file with the command:
      g++ -O2 -Wall -W -o unzcrash testsuite/unzcrash.cc
*/

#include <cstdio>
#include <cstdlib>
#include <stdint.h>
#include <signal.h>
#include <unistd.h>


void flip_bit( uint8_t buffer[], const int bit )
  {
  const int byteno = bit / 8;
  const int bitno  = bit % 8;
  const uint8_t mask = 1 << bitno;
//  std::fprintf( stderr, "(byte %d  bit %d  mask %d)", byteno, bitno, mask );
  buffer[byteno] ^= mask;
  }


int main( const int argc, const char * argv[] )
  {
  if( argc < 3 )
    {
    std::fprintf( stderr, "Usage: unzcrash \"lzip -tv\" filename.lz\n" );
    return 1;
    }

  FILE *f = std::fopen( argv[2], "rb" );
  if( !f )
    {
    std::fprintf( stderr, "Can't open file `%s' for reading\n", argv[2] );
    return 1;
    }

  const int buffer_size = 1 << 20;
  uint8_t buffer[buffer_size];
  const int size = std::fread( buffer, 1, buffer_size, f );
  if( size >= buffer_size )
    {
    std::fprintf( stderr, "input file `%s' too big.\n", argv[2] );
    return 1;
    }
  std::fclose( f );

  f = popen( argv[1], "w" );
  if( !f )
    {
    std::fprintf( stderr, "incorrect parameters or too many files.\n" );
    return 1;
    }
  const int wr = std::fwrite( buffer, 1, size, f );
  if( wr != size || pclose( f ) != 0 )
    {
    std::fprintf( stderr, "Could not run `%s' or other error.\n", argv[1] );
    return 1;
    }

  signal( SIGPIPE, SIG_IGN );

  for( int bit = 0; bit < size * 8; ++bit )
    {
    std::fprintf( stderr, "bit %d  ", bit );
    flip_bit( buffer, bit );
    f = popen( argv[1], "w" );
    if( !f )
      { std::fprintf( stderr, "Can't open pipe.\n" ); return 1; }
    std::fwrite( buffer, 1, size, f );
    pclose( f );
    flip_bit( buffer, bit );
    }

  return 0;
  }
