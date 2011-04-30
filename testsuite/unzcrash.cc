/*  Unzcrash - A test program written to test robustness to
               decompression of corrupted data.
    Inspired by unzcrash.c from Julian Seward's bzip2.
    Copyright (C) 2008, 2009, 2010, 2011 Antonio Diaz Diaz.

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <cerrno>
#include <climits>
#include <csignal>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>
#include <stdint.h>
#include <unistd.h>

#include "../arg_parser.h"

#if CHAR_BIT != 8
#error "Environments where CHAR_BIT != 8 are not supported."
#endif

#ifndef LLONG_MAX
#define LLONG_MAX  0x7FFFFFFFFFFFFFFFLL
#endif
#ifndef LLONG_MIN
#define LLONG_MIN  (-LLONG_MAX - 1LL)
#endif
#ifndef ULLONG_MAX
#define ULLONG_MAX 0xFFFFFFFFFFFFFFFFULL
#endif


namespace {

const char * const Program_name = "Unzcrash";
const char * const program_name = "unzcrash";
const char * const program_year = "2011";
const char * invocation_name = 0;

int verbosity = 0;


void show_help() throw()
  {
  std::printf( "%s - A test program written to test robustness to\n", Program_name );
  std::printf( "decompression of corrupted data.\n" );
  std::printf( "\nUsage: %s [options] \"lzip -tv\" filename.lz\n", invocation_name );
  std::printf( "\nThis program reads the specified file and then repeatedly decompresses\n" );
  std::printf( "it, increasing 256 times each byte of the compressed data, so as to test\n" );
  std::printf( "all possible one-byte errors. This should not cause any invalid memory\n" );
  std::printf( "accesses. If it does, please, report it as a bug.\n" );
  std::printf( "\nOptions:\n" );
  std::printf( "  -h, --help                 display this help and exit\n" );
  std::printf( "  -V, --version              output version information and exit\n" );
  std::printf( "  -b, --bits=<n>[,n]         test <n>-bit errors instead of full byte\n" );
  std::printf( "  -p, --position=<n>         first byte position to test\n" );
  std::printf( "  -q, --quiet                suppress all messages\n" );
  std::printf( "  -s, --size=<n>             number of byte positions to test\n" );
  std::printf( "  -v, --verbose              be verbose (a 2nd -v gives more)\n" );
  std::printf( "\nReport bugs to lzip-bug@nongnu.org\n" );
  std::printf( "Lzip home page: http://www.nongnu.org/lzip/lzip.html\n" );
  }


void show_version() throw()
  {
  std::printf( "%s %s\n", Program_name, PROGVERSION );
  std::printf( "Copyright (C) %s Antonio Diaz Diaz.\n", program_year );
  std::printf( "License GPLv3+: GNU GPL version 3 or later <http://gnu.org/licenses/gpl.html>\n" );
  std::printf( "This is free software: you are free to change and redistribute it.\n" );
  std::printf( "There is NO WARRANTY, to the extent permitted by law.\n" );
  }


void show_error( const char * const msg, const int errcode = 0,
                 const bool help = false ) throw()
  {
  if( verbosity >= 0 )
    {
    if( msg && msg[0] )
      {
      std::fprintf( stderr, "%s: %s", program_name, msg );
      if( errcode > 0 )
        std::fprintf( stderr, ": %s", std::strerror( errcode ) );
      std::fprintf( stderr, "\n" );
      }
    if( help && invocation_name && invocation_name[0] )
      std::fprintf( stderr, "Try `%s --help' for more information.\n",
                    invocation_name );
    }
  }


void internal_error( const char * const msg )
  {
  if( verbosity >= 0 )
    std::fprintf( stderr, "%s: internal error: %s.\n", program_name, msg );
  std::exit( 3 );
  }


long long getnum( const char * const ptr,
                  const long long llimit = LLONG_MIN + 1,
                  const long long ulimit = LLONG_MAX ) throw()
  {
  errno = 0;
  char *tail;
  long long result = strtoll( ptr, &tail, 0 );
  if( tail == ptr )
    {
    show_error( "Bad or missing numerical argument.", 0, true );
    std::exit( 1 );
    }

  if( !errno && tail[0] )
    {
    int factor = ( tail[1] == 'i' ) ? 1024 : 1000;
    int exponent = 0;
    bool bad_multiplier = false;
    switch( tail[0] )
      {
      case ' ': break;
      case 'Y': exponent = 8; break;
      case 'Z': exponent = 7; break;
      case 'E': exponent = 6; break;
      case 'P': exponent = 5; break;
      case 'T': exponent = 4; break;
      case 'G': exponent = 3; break;
      case 'M': exponent = 2; break;
      case 'K': if( factor == 1024 ) exponent = 1; else bad_multiplier = true;
                break;
      case 'k': if( factor == 1000 ) exponent = 1; else bad_multiplier = true;
                break;
      default : bad_multiplier = true;
      }
    if( bad_multiplier )
      {
      show_error( "Bad multiplier in numerical argument.", 0, true );
      std::exit( 1 );
      }
    for( int i = 0; i < exponent; ++i )
      {
      if( LLONG_MAX / factor >= llabs( result ) ) result *= factor;
      else { errno = ERANGE; break; }
      }
    }
  if( !errno && ( result < llimit || result > ulimit ) ) errno = ERANGE;
  if( errno )
    {
    show_error( "Numerical argument out of limits." );
    std::exit( 1 );
    }
  return result;
  }


class Bitset8			// 8 value bitset (1..8)
  {
  bool data[8];
  static bool valid_digit( const unsigned char ch ) throw()
    { return ( ch >= '1' && ch <= '8' ); }

public:
  Bitset8() throw() { for( int i = 0; i < 8; ++i ) data[i] = true; }

  bool includes( const int i ) const throw()
    { return ( i >= 1 && i <= 8 && data[i-1] ); }

  // Recognized formats: 1 1,2,3 1-4 1,3-5,8
  bool parse( const char * p ) throw()
    {
    for( int i = 0; i < 8; ++i ) data[i] = false;
    while( true )
      {
      const unsigned char ch1 = *p++;
      if( !valid_digit( ch1 ) ) break;
      if( *p != '-' ) data[ch1-'1'] = true;
      else
        {
        ++p;
        if( !valid_digit( *p ) || ch1 > *p ) break;
        for( int c = ch1; c <= *p; ++c ) data[c-'1'] = true;
        ++p;
        }
      if( *p == 0 ) return true;
      if( *p == ',' ) ++p; else break;
      }
    show_error( "Invalid value or range." );
    return false;
    }

  // number of n-bit errors per byte (n=0..8): 1 8 28 56 70 56 28 8 1
  void print() const throw()
    {
    std::fflush( stderr );
    int c = 0;
    for( int i = 0; i < 8; ++i ) if( data[i] ) ++c;
    if( c == 8 ) std::printf( "Testing full byte.\n" );
    else if( c == 0 ) std::printf( "Nothing to test.\n" );
    else
      {
      std::printf( "Testing " );
      for( int i = 0; i < 8; ++i )
        if( data[i] )
          {
          std::printf( "%d", i + 1 );
          if( --c ) std::printf( "," );
          }
      std::printf( " bit errors.\n" );
      }
    std::fflush( stdout );
    }
  };


int differing_bits( const uint8_t byte1, const uint8_t byte2 )
  {
  int count = 0;
  uint8_t dif = byte1 ^ byte2;
  while( dif )
    { count += ( dif & 1 ); dif >>= 1; }
  return count;
  }

} // end namespace


int main( const int argc, const char * const argv[] )
  {
  enum { buffer_size = 3 << 20 };
  Bitset8 bits;			// if Bitset8::parse not called test full byte
  int pos = 0;
  int max_size = buffer_size;
  invocation_name = argv[0];

  const Arg_parser::Option options[] =
    {
    { 'h', "help",     Arg_parser::no  },
    { 'b', "bits",     Arg_parser::yes },
    { 'p', "position", Arg_parser::yes },
    { 'q', "quiet",    Arg_parser::no  },
    { 's', "size",     Arg_parser::yes },
    { 'v', "verbose",  Arg_parser::no  },
    { 'V', "version",  Arg_parser::no  },
    {  0 ,  0,         Arg_parser::no  } };

  const Arg_parser parser( argc, argv, options );
  if( parser.error().size() )				// bad option
    { show_error( parser.error().c_str(), 0, true ); return 1; }

  int argind = 0;
  for( ; argind < parser.arguments(); ++argind )
    {
    const int code = parser.code( argind );
    if( !code ) break;					// no more options
    const char * const arg = parser.argument( argind ).c_str();
    switch( code )
      {
      case 'h': show_help(); return 0;
      case 'b': if( !bits.parse( arg ) ) return 1; break;
      case 'p': pos = getnum( arg, 0, buffer_size - 1 ); break;
      case 'q': verbosity = -1; break;
      case 's': max_size = getnum( arg, 1, buffer_size ); break;
      case 'v': if( verbosity < 4 ) ++verbosity; break;
      case 'V': show_version(); return 0;
      default : internal_error( "uncaught option" );
      }
    } // end process options

  if( argind + 2 != parser.arguments() )
    {
    if( verbosity >= 0 )
      std::fprintf( stderr, "Usage: %s \"lzip -tv\" filename.lz\n",
                    invocation_name );
    return 1;
    }

  FILE *f = std::fopen( parser.argument( argind + 1 ).c_str(), "rb" );
  if( !f )
    {
    if( verbosity >= 0 )
      std::fprintf( stderr, "Can't open file `%s' for reading\n",
                    parser.argument( argind + 1 ).c_str() );
    return 1;
    }

  uint8_t * const buffer = new uint8_t[buffer_size];
  const int size = std::fread( buffer, 1, buffer_size, f );
  if( size >= buffer_size )
    {
    if( verbosity >= 0 )
      std::fprintf( stderr, "input file `%s' is too big.\n",
                    parser.argument( argind + 1 ).c_str() );
    return 1;
    }
  std::fclose( f );

  f = popen( parser.argument( argind ).c_str(), "w" );
  if( !f )
    { show_error( "Can't open pipe", errno ); return 1; }
  const int wr = std::fwrite( buffer, 1, size, f );
  if( wr != size || pclose( f ) != 0 )
    {
    if( verbosity >= 0 )
      std::fprintf( stderr, "Could not run `%s' : %s.\n",
                    parser.argument( argind ).c_str(), std::strerror( errno ) );
    return 1;
    }

  std::signal( SIGPIPE, SIG_IGN );
  if( verbosity >= 1 ) bits.print();

  const int end = ( ( pos + max_size < size ) ? pos + max_size : size );
  for( int i = pos; i < end; ++i )
    {
    if( verbosity >= 0 )
      std::fprintf( stderr, "byte %d\n", i );
    const uint8_t byte = buffer[i];
    for( int j = 0; j < 255; ++j )
      {
      ++buffer[i];
      if( bits.includes( differing_bits( byte, buffer[i] ) ) )
        {
        f = popen( parser.argument( argind ).c_str(), "w" );
        if( !f )
          { show_error( "Can't open pipe", errno ); return 1; }
        std::fwrite( buffer, 1, size, f );
        pclose( f );
        }
      }
    buffer[i] = byte;
    }

  delete[] buffer;
  return 0;
  }
