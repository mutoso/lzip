/*  Lzip - A data compressor based on the LZMA algorithm
    Copyright (C) 2008, 2009 Antonio Diaz Diaz.

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

class Input_buffer
  {
  enum { buffer_size = 65536 };
  uint8_t * const buffer;
  int pos;
  int stream_pos;		// when reached, a new block must be read
  const int ides_;		// input file descriptor
  bool at_stream_end;

  bool read_block();

public:
  Input_buffer( const int ides )
    :
    buffer( new uint8_t[buffer_size] ),
    pos( 0 ),
    stream_pos( 0 ),
    ides_( ides ),
    at_stream_end( false ) {}

  ~Input_buffer() { delete[] buffer; }

  bool finished() const throw() { return at_stream_end && pos >= stream_pos; }

  uint8_t get_byte()
    {
    if( pos >= stream_pos && !read_block() ) return 0;
    return buffer[pos++];
    }
  };


class Range_decoder
  {
  mutable long long member_pos;
  uint32_t code;
  uint32_t range;
  Input_buffer & ibuf;

public:
  Range_decoder( const int header_size, Input_buffer & buf )
    :
    member_pos( header_size ),
    code( 0 ),
    range( 0xFFFFFFFF ),
    ibuf( buf )
    { for( int i = 0; i < 5; ++i ) code = (code << 8) | get_byte(); }

  bool finished() const throw() { return ibuf.finished(); }
  long long member_position() const throw() { return member_pos; }

  uint8_t get_byte() const
    {
    ++member_pos;
    return ibuf.get_byte();
    }

  void reload() throw()
    {
    code = 0;
    range = 0xFFFFFFFF;
    for( int i = 0; i < 5; ++i ) code = (code << 8) | get_byte();
    }

  void normalize()
    {
    if( range <= 0x00FFFFFF )
      { range <<= 8; code = (code << 8) | get_byte(); }
    }

  int decode( const int num_bits )
    {
    int symbol = 0;
    for( int i = num_bits; i > 0; --i )
      {
      symbol <<= 1;
      if( range <= 0x00FFFFFF )
        {
        range <<= 7; code = (code << 8) | get_byte();
        if( code >= range ) { code -= range; symbol |= 1; }
        }
      else
        {
        range >>= 1;
        if( code >= range ) { code -= range; symbol |= 1; }
        }
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
    int model = 1;
    for( int i = num_bits; i > 0; --i )
      model = ( model << 1 ) | decode_bit( bm[model] );
    return model - (1 << num_bits);
    }

  int decode_tree_reversed( Bit_model bm[], const int num_bits )
    {
    int model = 1;
    int symbol = 0;
    for( int i = 0; i < num_bits; ++i )
      {
      const int bit = decode_bit( bm[model] );
      model <<= 1;
      if( bit ) { model |= 1; symbol |= (1 << i); }
      }
    return symbol;
    }

  int decode_matched( Bit_model bm[], const int match_byte )
    {
    Bit_model *bm1 = bm + 0x100;
    int symbol = 1;
    for( int i = 1; i <= 8; ++i )
      {
      const int match_bit = ( match_byte << i ) & 0x100;
      const int bit = decode_bit( bm1[match_bit+symbol] );
      symbol = ( symbol << 1 ) | bit;
      if( ( match_bit && !bit ) || ( !match_bit && bit ) )
        {
        while( ++i <= 8 )
          symbol = ( symbol << 1 ) | decode_bit( bm[symbol] );
        break;
        }
      }
    return symbol & 0xFF;
    }
  };


class Len_decoder
  {
  Bit_model choice1;
  Bit_model choice2;
  Bit_model bm_low[pos_states][len_low_symbols];
  Bit_model bm_mid[pos_states][len_mid_symbols];
  Bit_model bm_high[len_high_symbols];

public:
  int decode( Range_decoder & range_decoder, const int pos_state )
    {
    if( range_decoder.decode_bit( choice1 ) == 0 )
      return range_decoder.decode_tree( bm_low[pos_state], len_low_bits );
    if( range_decoder.decode_bit( choice2 ) == 0 )
      return len_low_symbols +
             range_decoder.decode_tree( bm_mid[pos_state], len_mid_bits );
    return len_low_symbols + len_mid_symbols +
           range_decoder.decode_tree( bm_high, len_high_bits );
    }
  };


class Literal_decoder
  {
  Bit_model bm_literal[1<<literal_context_bits][0x300];

  int state( const int prev_byte ) const throw()
    { return ( prev_byte >> ( 8 - literal_context_bits ) ); }

public:
  uint8_t decode( Range_decoder & range_decoder, const uint8_t prev_byte )
    { return range_decoder.decode_tree( bm_literal[state(prev_byte)], 8 ); }

  uint8_t decode_matched( Range_decoder & range_decoder,
                          const uint8_t prev_byte, const uint8_t match_byte )
    { return range_decoder.decode_matched( bm_literal[state(prev_byte)], match_byte ); }
  };


class LZ_decoder
  {
  long long partial_data_pos;
  const int format_version;
  const int dictionary_size;
  const int buffer_size;
  uint8_t * const buffer;
  int pos;
  int stream_pos;		// first byte not yet written to file
  uint32_t crc_;
  const int odes_;		// output file descriptor

  Bit_model bm_match[State::states][pos_states];
  Bit_model bm_rep[State::states];
  Bit_model bm_rep0[State::states];
  Bit_model bm_rep1[State::states];
  Bit_model bm_rep2[State::states];
  Bit_model bm_len[State::states][pos_states];
  Bit_model bm_dis_slot[max_dis_states][1<<dis_slot_bits];
  Bit_model bm_dis[modeled_distances-end_dis_model];
  Bit_model bm_align[dis_align_size];

  Range_decoder range_decoder;
  Len_decoder len_decoder;
  Len_decoder rep_match_len_decoder;
  Literal_decoder literal_decoder;

  uint8_t get_byte( const int distance ) const throw()
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
      std::memcpy( buffer + pos, buffer + i, len );
      pos += len;
      }
    else for( ; len > 0 ; --len )
      {
      buffer[pos] = buffer[i];
      if( ++pos >= buffer_size ) flush_data();
      if( ++i >= buffer_size ) i = 0;
      }
    }

  void flush_data();
  bool verify_trailer( const Pretty_print & pp ) const;

public:
  LZ_decoder( const File_header & header, Input_buffer & ibuf, const int odes )
    :
    partial_data_pos( 0 ),
    format_version( header.version ),
    dictionary_size( header.dictionary_size() ),
    buffer_size( std::max( 65536, dictionary_size ) ),
    buffer( new uint8_t[buffer_size] ),
    pos( 0 ),
    stream_pos( 0 ),
    crc_( 0xFFFFFFFF ),
    odes_( odes ),
    range_decoder( sizeof header, ibuf ),
    literal_decoder()
    { buffer[buffer_size-1] = 0; }	// prev_byte of first_byte

  ~LZ_decoder() { delete[] buffer; }

  uint32_t crc() const throw() { return crc_ ^ 0xFFFFFFFF; }
  int decode_member( const Pretty_print & pp );

  long long member_position() const throw()
    { return range_decoder.member_position(); }
  long long data_position() const throw()
    { return partial_data_pos + pos; }
  };
