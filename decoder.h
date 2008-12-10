/*  Lzip - A data compressor based on the LZMA algorithm
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

class Input_buffer
  {
  const int buffer_size;
  uint8_t * const buffer;
  int pos;
  int size;		// when reached, a new block must be read
  const int _ides;
  bool _finished;

  bool read_block();

public:
  Input_buffer( const int ides )
    :
    buffer_size( 65536 ),
    buffer( new uint8_t[buffer_size] ),
    pos( 0 ),
    size( 0 ),
    _ides( ides ),
    _finished( false ) {}

  ~Input_buffer() { delete[] buffer; }

  bool finished() const throw() { return _finished; }

  uint8_t read_byte()
    {
    if( pos >= size && !read_block() ) return 0;
    return buffer[pos++];
    }
  };


class Range_decoder
  {
  long long file_pos;
  uint32_t code;
  uint32_t range;
  Input_buffer & ibuf;

public:
  Range_decoder( const int header_size, Input_buffer & buf )
    :
    file_pos( header_size ),
    code( 0 ), range( 0xFFFFFFFF ),
    ibuf( buf )
    { for( int i = 0; i < 5; ++i ) code = (code << 8) | read_byte(); }

  uint8_t read_byte()
    {
    ++file_pos;
    return ibuf.read_byte();
    }

  long long file_position() const throw() { return file_pos; }
  bool finished() const throw() { return ibuf.finished(); }

  int decode( const int num_bits = 1 )
    {
    int symbol = 0;
    for( int i = num_bits - 1; i >= 0; --i )
      {
      range >>= 1;
      symbol <<= 1;
      if( ( ( code - range ) & 0x80000000 ) == 0 )
        { code -= range; symbol |= 1; }
      if( range <= 0xFFFFFF )
        { range <<= 8; code = (code << 8) | read_byte(); }
      }
    return symbol;
    }

  int decode_bit( Bit_model & bm )
    {
    int symbol = 0;
    const uint32_t bound = ( range >> bit_model_total_bits ) * bm.probability;
    if( code < bound )
      {
      range = bound;
      bm.probability += (bit_model_total - bm.probability) >> bit_model_move_bits;
      }
    else
      {
      range -= bound;
      code -= bound;
      bm.probability -= bm.probability >> bit_model_move_bits;
      symbol = 1;
      }
    if( range <= 0xFFFFFF )
      { range <<= 8; code = (code << 8) | read_byte(); }
    return symbol;
    }

  int decode_tree( Bit_model bm[], const int num_bits )
    {
    int model = 1;
    for( int i = num_bits; i > 0; --i )
      model = ( model << 1 ) | decode_bit( bm[model-1] );
    return model - (1 << num_bits);
    }

  int decode_tree_reversed( Bit_model bm[], const int num_bits )
    {
    int model = 1;
    int symbol = 0;
    for( int i = 0; i < num_bits; ++i )
      {
      model = ( model << 1 ) | decode_bit( bm[model-1] );
      if( model & 1 ) symbol |= ( 1 << i );
      }
    return symbol;
    }

  int decode_matched( Bit_model bm[], const int match_byte )
    {
    int symbol = 1;
    for( int i = 7; i >= 0; --i )
      {
      int match_bit = ( match_byte >> i ) & 1;
      int bit = decode_bit( bm[(match_bit<<8)+symbol+0xFF] );
      symbol = ( symbol << 1 ) | bit;
      if( match_bit != bit ) break;
      }
    while( symbol < 0x100 )
      symbol = ( symbol << 1 ) | decode_bit( bm[symbol-1] );
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
  typedef Bit_model Bm_array[0x300];
  Bm_array * const bm_literal;

  int state( const int prev_byte ) const throw()
    { return ( prev_byte >> ( 8 - literal_context_bits ) ); }

public:
  Literal_decoder()
    : bm_literal( new Bm_array[1<<literal_context_bits] ) {}

  ~Literal_decoder() { delete[] bm_literal; }

  uint8_t decode( Range_decoder & range_decoder, const int prev_byte )
    { return range_decoder.decode_tree( bm_literal[state(prev_byte)], 8 ); }

  uint8_t decode_matched( Range_decoder & range_decoder,
                          const int prev_byte, const int match_byte )
    { return range_decoder.decode_matched( bm_literal[state(prev_byte)], match_byte ); }
  };


class LZ_decoder
  {
  long long partial_file_pos;
  const int buffer_size;
  uint8_t * const buffer;
  int pos;
  int stream_pos;
  const int _odes;
  uint32_t _crc;

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
    int newpos = pos - distance - 1;
    if( newpos < 0 ) newpos += buffer_size;
    return buffer[newpos];
    }

  void put_byte( const uint8_t b )
    {
    update_crc( _crc, b );
    buffer[pos] = b;
    if( ++pos >= buffer_size ) flush();
    }

  bool copy_block( const int distance, int len )
    {
    if( distance < 0 || distance >= buffer_size ||
        len <= 0 || len > max_match_len ) return false;
    int newpos = pos - distance - 1;
    if( newpos < 0 ) newpos += buffer_size;
    for( ; len > 0 ; --len )
      {
      update_crc( _crc, buffer[newpos] );
      buffer[pos] = buffer[newpos];
      if( ++pos >= buffer_size ) flush();
      if( ++newpos >= buffer_size ) newpos = 0;
      }
    return true;
    }

  void flush();
  bool verify_trailer( const Pretty_print & pp );

public:
  LZ_decoder( const File_header & header, Input_buffer & ibuf, const int odes )
    :
    partial_file_pos( 0 ),
    buffer_size( 1 << header.dictionary_bits ),
    buffer( new uint8_t[buffer_size] ),
    pos( 0 ),
    stream_pos( 0 ),
    _odes( odes ),
    _crc( 0xFFFFFFFF ),
    range_decoder( sizeof header, ibuf ),
    literal_decoder() {}

  ~LZ_decoder() { delete[] buffer; }

  uint32_t crc() const throw() { return _crc ^ 0xFFFFFFFF; }
  int decode( const Pretty_print & pp );

  long long input_file_position() const throw()
    { return range_decoder.file_position(); }
  long long output_file_position() const throw()
    { return partial_file_pos + pos; }
  };
