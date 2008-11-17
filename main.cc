/*  Lzip - A LZMA file compressor
    Copyright (C) 2008 Antonio Diaz Diaz.

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
/*
    Return values: 0 for a normal exit, 1 for environmental problems
    (file not found, invalid flags, I/O errors, etc), 2 to indicate a
    corrupt or invalid input file, 3 for an internal consistency error
    (eg, bug) which caused lzip to panic.
*/

#define _FILE_OFFSET_BITS 64

#include <algorithm>
#include <cerrno>
#include <climits>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>
#include <fcntl.h>
#include <stdint.h>
#include <signal.h>
#include <unistd.h>
#include <utime.h>
#include <sys/stat.h>
#include <sys/time.h>

#include "arg_parser.h"
#include "lzip.h"
#include "decoder.h"
#include "encoder.h"

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

const char * invocation_name = 0;
const char * const Program_name    = "Lzip";
const char * const program_name    = "lzip";
const char * const program_year    = "2008";

struct { const char * from; const char * to; } const known_extensions[] = {
  { ".lz",  ""     },
  { ".tlz", ".tar" },
  { 0,      0      } };

struct lzma_options
  {
  int dictionary_bits;		// 12..30
  int match_len_limit;		// 5..273
  };

enum Mode { m_compress = 0, m_decompress, m_test };

std::string input_filename;
std::string output_filename;
int outhandle = -1;
bool delete_output_on_interrupt = false;


void show_help() throw()
  {
  std::printf( "%s - A LZMA file compressor.\n", Program_name );
  std::printf( "\nUsage: %s [options] [files]\n", invocation_name );
  std::printf( "Options:\n" );
  std::printf( "  -h, --help                 display this help and exit\n" );
  std::printf( "  -V, --version              output version information and exit\n" );
  std::printf( "  -c, --stdout               send output to standard output\n" );
  std::printf( "  -d, --decompress           force decompression\n" );
  std::printf( "  -f, --force                overwrite existing output files\n" );
  std::printf( "  -k, --keep                 keep (don't delete) input files\n" );
  std::printf( "  -m, --match-length=<n>     set match length limit in bytes [64]\n" );
  std::printf( "  -q, --quiet                suppress all messages\n" );
  std::printf( "  -s, --dictionary-size=<n>  set dictionary size in bytes [8MiB]\n" );
  std::printf( "  -t, --test                 test compressed file integrity\n" );
  std::printf( "  -v, --verbose              be verbose (a 2nd -v gives more)\n" );
  std::printf( "  -z, --compress             force compression\n" );
  std::printf( "  -1 .. -9                   set compression level [default 6]\n" );
  std::printf( "      --fast                 alias for -1\n" );
  std::printf( "      --best                 alias for -9\n" );
  std::printf( "If no file names are given, lzip compresses or decompresses\n" );
  std::printf( "from standard input to standard output.\n" );
  std::printf( "\nReport bugs to lzip-bug@nongnu.org\n");
  }


void show_version() throw()
  {
  std::printf( "%s %s\n", Program_name, PROGVERSION );
  std::printf( "Copyright (C) %s Antonio Diaz Diaz.\n", program_year );
  std::printf( "License GPLv3+: GNU GPL version 3 or later <http://gnu.org/licenses/gpl.html>\n" );
  std::printf( "This is free software: you are free to change and redistribute it.\n" );
  std::printf( "There is NO WARRANTY, to the extent permitted by law.\n" );
  }


const char * format_num( long long num, long long max = 1023,
                         const int set_prefix = 0 ) throw()
  {
  const char * const si_prefix[8] =
    { "k", "M", "G", "T", "P", "E", "Z", "Y" };
  const char * const binary_prefix[8] =
    { "Ki", "Mi", "Gi", "Ti", "Pi", "Ei", "Zi", "Yi" };
  static bool si = false;
  static char buf[16];

  if( set_prefix ) si = ( set_prefix > 0 );
  const int factor = ( si ) ? 1000 : 1024;
  const char * const *prefix = ( si ) ? si_prefix : binary_prefix;
  const char *p = "";
  max = std::max( 999LL, std::min( 999999LL, max ) );

  for( int i = 0; i < 8 && llabs( num ) > llabs( max ); ++i )
    { num /= factor; p = prefix[i]; }
  snprintf( buf, sizeof buf, "%lld %s", num, p );
  return buf;
  }


long long getnum( const char * ptr, const int bs,
                  const long long min = LLONG_MIN + 1,
                  const long long max = LLONG_MAX ) throw()
  {
  errno = 0;
  char *tail;
  long long result = strtoll( ptr, &tail, 0 );
  if( tail == ptr )
    {
    show_error( "bad or missing numerical argument", 0, true );
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
      case 'b': if( bs > 0 ) { factor = bs; exponent = 1; }
                else bad_multiplier = true;
                break;
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
      default: bad_multiplier = true;
      }
    if( bad_multiplier )
      {
      show_error( "bad multiplier in numerical argument", 0, true );
      std::exit( 1 );
      }
    for( int i = 0; i < exponent; ++i )
      {
      if( LLONG_MAX / factor >= llabs( result ) ) result *= factor;
      else { errno = ERANGE; break; }
      }
    }
  if( !errno && ( result < min || result > max ) ) errno = ERANGE;
  if( errno )
    {
    show_error( "numerical argument out of limits" );
    std::exit( 1 );
    }
  return result;
  }


int real_bits( const int value ) throw()
  {
  int bits = 0;
  for( int i = 1, mask = 1; mask > 0; ++i, mask <<= 1 )
    if( value & mask ) bits = i;
  return bits;
  }


int get_dict_bits( const char * arg ) throw()
  {
  char *tail;
  int bits = std::strtol( arg, &tail, 0 );
  if( bits >= min_dictionary_bits && bits <= max_dictionary_bits && *tail == 0 )
    return bits;
  const int min_size = 1 << min_dictionary_bits;
  const int max_size = 1 << max_dictionary_bits;
  int size = getnum( arg, 0, ( min_size / 2 ) + 1, max_size );
  bits = real_bits( size - 1 );
  if( size > ( 1 << bits ) ) ++bits;
  return bits;
  }


int compress( const int ides, const int odes, lzma_options encoder_options,
              const Pretty_print & pp ) throw()
  {
  File_header header;
  header.set_magic();
  header.dictionary_bits = encoder_options.dictionary_bits;
  const int rd = writeblock( odes, (char *)&header, sizeof header );
  if( rd != sizeof header )
    { pp(); show_error( "error writing file header", errno ); return 1; }

  try {
    if( verbosity >= 1 ) pp();
    LZ_encoder encoder( header, ides, odes, encoder_options.match_len_limit );

    if( !encoder.encode() ) { pp( "encoder error" ); return 2; }

    if( verbosity >= 1 )
      {
      long long in_size = encoder.input_file_position();
      long long out_size = encoder.output_file_position();

      if( in_size <= 0 || out_size <= 0 )
        std::fprintf( stderr, "no data compressed.\n" );
      else
        std::fprintf( stderr, "%6.3f:1, %6.3f bits/byte, "
                              "%5.2f%% saved, %lld in, %lld out.\n",
                      (double)in_size / out_size,
                      ( 8.0 * out_size ) / in_size,
                      100.0 * ( 1.0 - ( (double)out_size / in_size ) ),
                      in_size, out_size );
      }
    }
  catch( std::bad_alloc )
    {
    pp( "not enough memory. Try a smaller dictionary size" );
    return 1;
    }
  catch( Error e ) { pp(); show_error( e.s, errno ); return 1; }
  return 0;
  }


int decompress( const int ides, const int odes, const Pretty_print & pp,
                const bool testing ) throw()
  {
  Input_buffer ibuf( ides );
  for( bool first_pass = true; ; first_pass = false, pp.reset() )
    {
    File_header header;
    for( unsigned int i = 0; i < sizeof header; ++i )
      ((uint8_t *)&header)[i] = ibuf.read_byte();
    if( ibuf.finished() )
      {
      if( first_pass ) { pp( "error reading file header" ); return 1; }
      else break;
      }
    if( !header.verify_magic() )
      {
      if( !first_pass ) break;
      if( verbosity >= 0 )
        { pp();
          std::fprintf( stderr, "bad magic number (file not created by %s).\n",
                        program_name ); }
      return 2;
      }
    if( !header.verify_version() )
      {
      if( verbosity >= 0 )
        { pp();
          std::fprintf( stderr, "file format not supported, newer %s needed.\n",
                        program_name ); }
      return 2;
      }
    if( header.dictionary_bits < min_dictionary_bits ||
        header.dictionary_bits > max_dictionary_bits )
      { pp( "invalid value in file header" ); return 2; }

    try {
      if( verbosity >= 1 )
        {
        pp();
        if( verbosity >= 2 )
          std::fprintf( stderr, "version %d, dictionary size %6sB.  ",
                        header.version,
                        format_num( 1 << header.dictionary_bits ) );
        }
      LZ_decoder decoder( header, ibuf, odes );

      const int result = decoder.decode( pp );
      if( result != 0 )
        {
        if( verbosity >= 0 && result <= 2 )
          {
          pp();
          if( result == 1 )
            std::fprintf( stderr, "decoder error at pos %lld\n",
                          decoder.input_file_position() );
          if( result == 2 )
            std::fprintf( stderr, "file ends unexpectedly at pos %lld\n",
                          decoder.input_file_position() );
          }
        return 2;
        }
      if( verbosity >= 1 )
        { if( testing ) std::fprintf( stderr, "ok\n" );
          else std::fprintf( stderr, "done\n" ); }
      }
    catch( std::bad_alloc )
      {
      pp( "not enough memory. Find a machine with more memory" );
      return 1;
      }
    catch( Error e ) { pp(); show_error( e.s, errno ); return 1; }
    }
  return 0;
  }


int extension_index() throw()
  {
  for( int i = 0; known_extensions[i].from; ++i )
    {
    const std::string ext( known_extensions[i].from );
    if( input_filename.size() > ext.size() &&
        input_filename.compare( input_filename.size() - ext.size(),
                                ext.size(), ext ) == 0 )
      return i;
    }
  return -1;
  }


std::string replace_extension( const int i ) throw()
  {
  if( i >= 0 )
    {
    const std::string from = known_extensions[i].from;
    if( input_filename.size() > from.size() )
      {
      const int len = input_filename.size() - from.size();
      std::string out( input_filename, 0, len );
      out += known_extensions[i].to;
      return out;
      }
    }
  std::string out( input_filename ); out += ".out";
  if( verbosity >= 0 )
    std::fprintf( stderr, "%s: can't guess original name for `%s' -- using `%s'.\n",
                  program_name, input_filename.c_str(), out.c_str() );
  return out;
  }


int open_instream( struct stat * in_statsp, const Mode program_mode,
                   const int eindex, const bool force ) throw()
  {
  if( program_mode == m_compress && !force && eindex >= 0 )
    {
    if( verbosity >= 0 )
      std::fprintf( stderr, "%s: input file `%s' already has `%s' suffix.\n",
                    program_name, input_filename.c_str(),
                    known_extensions[eindex].from );
    return -1;
    }

  int ides = open( input_filename.c_str(), O_RDONLY );
  if( ides < 0 )
    {
    if( verbosity >= 0 )
      std::fprintf( stderr, "%s: Can't open input file `%s': %s.\n",
                    program_name, input_filename.c_str(), strerror( errno ) );
    }
  else
    {
    const int i = fstat( ides, in_statsp );
    if( i < 0 || !S_ISREG( in_statsp->st_mode ) )
      {
      if( verbosity >= 0 )
        std::fprintf( stderr, "%s: input file `%s' is not a regular file.\n",
                      program_name, input_filename.c_str() );
      close( ides );
      ides = -1;
      }
    }
  return ides;
  }


int open_outstream( const Mode program_mode, const int eindex, const bool force ) throw()
  {
  int odes;

  if( program_mode == m_test )
    output_filename = "/dev/null";
  else if( program_mode == m_decompress )
    output_filename = replace_extension( eindex );
  else
    output_filename = input_filename + known_extensions[0].from;

  if( force || ( program_mode == m_test ) )
    odes = open( output_filename.c_str(), O_CREAT | O_TRUNC | O_WRONLY, S_IRUSR | S_IWUSR );
  else odes = open( output_filename.c_str(), O_CREAT | O_EXCL | O_WRONLY, S_IRUSR | S_IWUSR );
  if( odes < 0 && verbosity >= 0 )
    {
    if( errno == EEXIST )
      std::fprintf( stderr, "%s: Output file %s already exists, skipping.\n",
                    program_name, output_filename.c_str() );
    else
      std::fprintf( stderr, "%s: Can't create output file `%s': %s.\n",
                    program_name, output_filename.c_str(), strerror( errno ) );
    }
  return odes;
  }


bool check_tty( const Mode program_mode, const int inhandle ) throw()
  {
  if( program_mode == m_compress && isatty( outhandle ) )
    {
    show_error( "I won't write compressed data to a terminal.", 0, true );
    return false;
    }
  if( ( program_mode == m_decompress || program_mode == m_test ) &&
      isatty( inhandle ) )
    {
    show_error( "I won't read compressed data from a terminal.", 0, true );
    return false;
    }
  return true;
  }


void cleanup_and_fail( const int retval ) throw()
  {
  if( delete_output_on_interrupt )
    {
    if( verbosity >= 0 )
      std::fprintf( stderr, "%s: Deleting output file `%s', if it exists.\n",
               program_name, output_filename.c_str() );
    if( outhandle >= 0 ) { close( outhandle ); outhandle = -1; }
    if( std::remove( output_filename.c_str() ) != 0 )
      show_error( "WARNING: deletion of output file (apparently) failed." );
    }
  std::exit( retval );
  }


void signal_handler( const int ) throw()
  {
  show_error( "Control-C or similar caught, quitting." );
  cleanup_and_fail( 0 );
  }


void set_signals() throw()
  {
  signal( SIGTERM, signal_handler );
  signal( SIGHUP, signal_handler );
  signal( SIGINT, signal_handler );
  }

} // end namespace


int verbosity = 0;


void Pretty_print::operator()( const char * const msg ) const throw()
  {
  if( verbosity >= 0 )
    {
    if( first_post )
      {
      first_post = false;
      std::fprintf( stderr, "  %s: ", _name.c_str() );
      for( unsigned int i = 0; i < longest_name - _name.size(); ++i )
        std::fprintf( stderr, " " );
      if( !msg ) std::fflush( stderr );
      }
    if( msg ) std::fprintf( stderr, "%s.\n", msg );
    }
  }


void show_error( const char * msg, const int errcode, const bool help ) throw()
  {
  if( verbosity >= 0 )
    {
    if( msg && msg[0] != 0 )
      {
      std::fprintf( stderr, "%s: %s", program_name, msg );
      if( errcode > 0 ) std::fprintf( stderr, ": %s", strerror( errcode ) );
      std::fprintf( stderr, "\n" );
      }
    if( help && invocation_name && invocation_name[0] != 0 )
      std::fprintf( stderr, "Try `%s --help' for more information.\n", invocation_name );
    }
  }


void internal_error( const char * msg ) throw()
  {
  std::string s( "internal error: " ); s += msg;
  show_error( s.c_str() );
  std::exit( 3 );
  }


// Returns the number of bytes really read.
// If (returned value < size) and (errno == 0), means EOF was reached.
//
int readblock( const int fd, char * buf, const int size ) throw()
  {
  int rest = size;
  errno = 0;
  while( rest > 0 )
    {
    errno = 0;
    const int n = read( fd, buf + size - rest, rest );
    if( n > 0 ) rest -= n;
    else if( n == 0 ) break;
    else if( errno != EINTR && errno != EAGAIN ) break;
    }
  return ( rest > 0 ) ? size - rest : size;
  }


// Returns the number of bytes really written.
// If (returned value < size), it is always an error.
//
int writeblock( const int fd, const char * buf, const int size ) throw()
  {
  int rest = size;
  errno = 0;
  while( rest > 0 )
    {
    errno = 0;
    const int n = write( fd, buf + size - rest, rest );
    if( n > 0 ) rest -= n;
    else if( errno && errno != EINTR && errno != EAGAIN ) break;
    }
  return ( rest > 0 ) ? size - rest : size;
  }


int main( const int argc, const char * argv[] ) throw()
  {
  // Mapping from gzip/bzip2 style 1..9 compression modes
  // to the corresponding LZMA compression modes.
  const lzma_options option_mapping[] =
    {
    { 22,  10 },		// -1
    { 22,  12 },		// -2
    { 22,  16 },		// -3
    { 22,  32 },		// -4
    { 22,  64 },		// -5
    { 23,  64 },		// -6
    { 24,  64 },		// -7
    { 24, 128 },		// -8
    { 25, 273 } };		// -9
  lzma_options encoder_options = option_mapping[5];	// default = "-6"
  Mode program_mode = m_compress;
  int inhandle = -1;
  bool force = false;
  bool keep_input_files = false;
  bool to_stdout = false;
  std::vector< std::string > filenames;
  invocation_name = argv[0];

  const Arg_parser::Option options[] =
    {
    { '1', "fast",            Arg_parser::no  },
    { '2',  0,                Arg_parser::no  },
    { '3',  0,                Arg_parser::no  },
    { '4',  0,                Arg_parser::no  },
    { '5',  0,                Arg_parser::no  },
    { '6',  0,                Arg_parser::no  },
    { '7',  0,                Arg_parser::no  },
    { '8',  0,                Arg_parser::no  },
    { '9', "best",            Arg_parser::no  },
    { 'c', "stdout",          Arg_parser::no  },
    { 'd', "decompress",      Arg_parser::no  },
    { 'f', "force",           Arg_parser::no  },
    { 'h', "help",            Arg_parser::no  },
    { 'k', "keep",            Arg_parser::no  },
    { 'm', "match-length",    Arg_parser::yes },
    { 'q', "quiet",           Arg_parser::no  },
    { 's', "dictionary-size", Arg_parser::yes },
    { 't', "test",            Arg_parser::no  },
    { 'v', "verbose",         Arg_parser::no  },
    { 'V', "version",         Arg_parser::no  },
    { 'z', "compress",        Arg_parser::no  },
    {  0 ,  0,                Arg_parser::no  } };

  Arg_parser parser( argc, argv, options );
  if( parser.error().size() )				// bad option
    { show_error( parser.error().c_str(), 0, true ); return 1; }

  int argind = 0;
  for( ; argind < parser.arguments(); ++argind )
    {
    const int code = parser.code( argind );
    if( !code ) break;					// no more options
    const char * arg = parser.argument( argind ).c_str();
    switch( code )
      {
      case '1': case '2': case '3':
      case '4': case '5': case '6':
      case '7': case '8': case '9':
                encoder_options = option_mapping[code-'1']; break;
      case 'c': to_stdout = true; break;
      case 'd': program_mode = m_decompress; break;
      case 'f': force = true; break;
      case 'h': show_help(); return 0;
      case 'k': keep_input_files = true; break;
      case 'm': encoder_options.match_len_limit =
                getnum( arg, 0, 5, max_match_len ); break;
      case 'q': verbosity = -1; break;
      case 's': encoder_options.dictionary_bits = get_dict_bits( arg );
                break;
      case 't': program_mode = m_test; break;
      case 'v': if( verbosity < 4 ) ++verbosity; break;
      case 'V': show_version(); return 0;
      case 'z': program_mode = m_compress; break;
      default : internal_error( "uncaught option" );
      }
    }

  bool filenames_given = false;
  for( ; argind < parser.arguments(); ++argind )
    {
    if( parser.argument( argind ) != "-" ) filenames_given = true;
    filenames.push_back( parser.argument( argind ) );
    }

  if( filenames.empty() ) filenames.push_back("-");
  if( filenames_given ) set_signals();

  Pretty_print pp( filenames );
  if( program_mode == m_test )
    {
    outhandle = open_outstream( program_mode, 0, force );
    if( outhandle < 0 ) return 1;
    }

  int retval = 0;
  for( unsigned int i = 0; i < filenames.size(); ++i )
    {
    struct stat in_stats;
    output_filename.clear();

    if( !filenames[i].size() || filenames[i] == "-" )
      {
      input_filename.clear();
      inhandle = STDIN_FILENO;
      if( program_mode != m_test ) outhandle = STDOUT_FILENO;
      }
    else
      {
      input_filename = filenames[i];
      const int eindex = extension_index();
      inhandle = open_instream( &in_stats, program_mode, eindex, force );
      if( inhandle < 0 ) continue;
      if( program_mode != m_test )
        {
        if( to_stdout ) outhandle = STDOUT_FILENO;
        else
          {
          outhandle = open_outstream( program_mode, eindex, force );
          if( outhandle < 0 ) return 1;
          }
        }
      }

    if( !check_tty( program_mode, inhandle ) ) return 1;

    if( input_filename.size() && !to_stdout && program_mode != m_test )
      delete_output_on_interrupt = true;
    pp.set_name( input_filename );
    int tmp = 0;
    if( program_mode == m_compress )
      tmp = compress( inhandle, outhandle, encoder_options, pp );
    else
      tmp = decompress( inhandle, outhandle, pp, program_mode == m_test );
    if( tmp > retval ) retval = tmp;
    if( retval && program_mode != m_test ) cleanup_and_fail( retval );

    // Set permissions, owner and times.
    if( input_filename.size() && !to_stdout && program_mode != m_test )
      {
      tmp = 0;
      if( fchmod( outhandle, in_stats.st_mode ) != 0 ) tmp = 1;
      if( !tmp ) fchown( outhandle, in_stats.st_uid, in_stats.st_gid );
      if( close( outhandle ) != 0 ) cleanup_and_fail( 1 );
      delete_output_on_interrupt = false;
      if( !tmp )
	{
        struct utimbuf t;
        t.actime = in_stats.st_atime;
        t.modtime = in_stats.st_mtime;
        tmp = utime( output_filename.c_str(), &t );
        }
      if( tmp )
	{
        if( tmp > retval ) retval = tmp;
        show_error( "I can't change output file attributes." );
        cleanup_and_fail( retval );
	}
      }
    if( input_filename.size() )
      {
      close( inhandle ); inhandle = -1;
      if( !keep_input_files && !to_stdout && program_mode != m_test )
        std::remove( input_filename.c_str() );
      }
    }
  return retval;
  }
