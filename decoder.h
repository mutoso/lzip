/*  Lzip - LZMA lossless data compressor
    Copyright (C) 2008, 2009, 2010, 2011, 2012, 2013 Antonio Diaz Diaz.

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

class Range_decoder
  {
  enum { buffer_size = 16384 };
  unsigned long long partial_member_pos;
  uint8_t * const buffer;	// input buffer
  int pos;			// current pos in buffer
  int stream_pos;		// when reached, a new block must be read
  uint32_t code;
  uint32_t range;
  const int infd;		// input file descriptor
  bool at_stream_end;

  bool read_block();

  Range_decoder( const Range_decoder & );	// declared as private
  void operator=( const Range_decoder & );	// declared as private

public:
  explicit Range_decoder( const int ifd )
    :
    partial_member_pos( 0 ),
    buffer( new uint8_t[buffer_size] ),
    pos( 0 ),
    stream_pos( 0 ),
    code( 0 ),
    range( 0xFFFFFFFFU ),
    infd( ifd ),
    at_stream_end( false )
    {}

  ~Range_decoder() { delete[] buffer; }

  bool code_is_zero() const { return ( code == 0 ); }
  bool finished() { return pos >= stream_pos && !read_block(); }
  unsigned long long member_position() const { return partial_member_pos + pos; }
  void reset_member_position() { partial_member_pos = -pos; }

  uint8_t get_byte()
    {
    if( finished() ) return 0xAA;		// make code != 0
    return buffer[pos++];
    }

  int read_data( uint8_t * const outbuf, const int size )
    {
    int rest = size;
    while( rest > 0 && !finished() )
      {
      const int rd = std::min( rest, stream_pos - pos );
      std::memcpy( outbuf + size - rest, buffer + pos, rd );
      pos += rd;
      rest -= rd;
      }
    return size - rest;
    }

  void load()
    {
    code = 0;
    for( int i = 0; i < 5; ++i ) code = (code << 8) | get_byte();
    range = 0xFFFFFFFFU;
    }

  void normalize()
    {
    if( range <= 0x00FFFFFFU )
      { range <<= 8; code = (code << 8) | get_byte(); }
    }

  int decode( const int num_bits )
    {
    int symbol = 0;
    for( int i = num_bits; i > 0; --i )
      {
      normalize();
      range >>= 1;
//      symbol <<= 1;
//      if( code >= range ) { code -= range; symbol |= 1; }
      const uint32_t mask = 0U - (code < range);
      code -= range;
      code += range & mask;
      symbol = (symbol << 1) + (mask + 1);
      }
    return symbol;
    }

  int decode_bit( Bit_model & bm )
    {
    normalize();
    const uint32_t bound = ( range >> bit_model_total_bits ) * bm.probability;
    if( code < bound )
      {
      range = bound;
      bm.probability += (bit_model_total - bm.probability) >> bit_model_move_bits;
      return 0;
      }
    else
      {
      range -= bound;
      code -= bound;
      bm.probability -= bm.probability >> bit_model_move_bits;
      return 1;
      }
    }

  int decode_tree( Bit_model bm[], const int num_bits )
    {
    int symbol = 1;
    for( int i = num_bits; i > 0; --i )
      symbol = ( symbol << 1 ) | decode_bit( bm[symbol] );
    return symbol - (1 << num_bits);
    }

  int decode_tree6( Bit_model bm[] )
    {
    int symbol = 1;
    symbol = ( symbol << 1 ) | decode_bit( bm[symbol] );
    symbol = ( symbol << 1 ) | decode_bit( bm[symbol] );
    symbol = ( symbol << 1 ) | decode_bit( bm[symbol] );
    symbol = ( symbol << 1 ) | decode_bit( bm[symbol] );
    symbol = ( symbol << 1 ) | decode_bit( bm[symbol] );
    symbol = ( symbol << 1 ) | decode_bit( bm[symbol] );
    return symbol - (1 << 6);
    }

  int decode_tree_reversed( Bit_model bm[], const int num_bits )
    {
    int model = 1;
    int symbol = 0;
    for( int i = 0; i < num_bits; ++i )
      {
      const bool bit = decode_bit( bm[model] );
      model <<= 1;
      if( bit ) { ++model; symbol |= (1 << i); }
      }
    return symbol;
    }

  int decode_tree_reversed4( Bit_model bm[] )
    {
    int model = 1;
    int symbol = 0;
    int bit = decode_bit( bm[model] );
    model = (model << 1) + bit; symbol |= bit;
    bit = decode_bit( bm[model] );
    model = (model << 1) + bit; symbol |= (bit << 1);
    bit = decode_bit( bm[model] );
    model = (model << 1) + bit; symbol |= (bit << 2);
    if( decode_bit( bm[model] ) ) symbol |= 8;
    return symbol;
    }

  int decode_matched( Bit_model bm[], int match_byte )
    {
    Bit_model * const bm1 = bm + 0x100;
    int symbol = 1;
    for( int i = 7; i >= 0; --i )
      {
      match_byte <<= 1;
      const int match_bit = match_byte & 0x100;
      const int bit = decode_bit( bm1[match_bit+symbol] );
      symbol = ( symbol << 1 ) | bit;
      if( match_bit != bit << 8 )
        {
        while( symbol < 0x100 )
          symbol = ( symbol << 1 ) | decode_bit( bm[symbol] );
        break;
        }
      }
    return symbol - 0x100;
    }

  int decode_len( Len_model & lm, const int pos_state )
    {
    if( decode_bit( lm.choice1 ) == 0 )
      return decode_tree( lm.bm_low[pos_state], len_low_bits );
    if( decode_bit( lm.choice2 ) == 0 )
      return len_low_symbols +
             decode_tree( lm.bm_mid[pos_state], len_mid_bits );
    return len_low_symbols + len_mid_symbols +
           decode_tree( lm.bm_high, len_high_bits );
    }
  };


class LZ_decoder
  {
  unsigned long long partial_data_pos;
  Range_decoder & rdec;
  const int dictionary_size;
  const int buffer_size;
  uint8_t * const buffer;	// output buffer
  int pos;			// current pos in buffer
  int stream_pos;		// first byte not yet written to file
  uint32_t crc_;
  const int outfd;		// output file descriptor
  const int member_version;

  void flush_data();
  bool verify_trailer( const Pretty_print & pp ) const;

  uint8_t get_prev_byte() const
    {
    const int i = ( ( pos > 0 ) ? pos : buffer_size ) - 1;
    return buffer[i];
    }

  uint8_t get_byte( const int distance ) const
    {
    int i = pos - distance - 1;
    if( i < 0 ) i += buffer_size;
    return buffer[i];
    }

  void put_byte( const uint8_t b )
    {
    buffer[pos] = b;
    if( ++pos >= buffer_size ) flush_data();
    }

  void copy_block( const int distance, int len )
    {
    int i = pos - distance - 1;
    if( i < 0 ) i += buffer_size;
    if( len < buffer_size - std::max( pos, i ) && len <= std::abs( pos - i ) )
      {
      std::memcpy( buffer + pos, buffer + i, len );	// no wrap, no overlap
      pos += len;
      }
    else for( ; len > 0; --len )
      {
      buffer[pos] = buffer[i];
      if( ++pos >= buffer_size ) flush_data();
      if( ++i >= buffer_size ) i = 0;
      }
    }

  LZ_decoder( const LZ_decoder & );		// declared as private
  void operator=( const LZ_decoder & );		// declared as private

public:
  LZ_decoder( const File_header & header, Range_decoder & rde, const int ofd )
    :
    partial_data_pos( 0 ),
    rdec( rde ),
    dictionary_size( header.dictionary_size() ),
    buffer_size( std::max( 65536, dictionary_size ) ),
    buffer( new uint8_t[buffer_size] ),
    pos( 0 ),
    stream_pos( 0 ),
    crc_( 0xFFFFFFFFU ),
    outfd( ofd ),
    member_version( header.version() )
    { buffer[buffer_size-1] = 0; }	// prev_byte of first_byte

  ~LZ_decoder() { delete[] buffer; }

  unsigned crc() const { return crc_ ^ 0xFFFFFFFFU; }
  unsigned long long data_position() const { return partial_data_pos + pos; }

  int decode_member( const Pretty_print & pp );
  };
