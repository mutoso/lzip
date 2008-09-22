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

const int max_num_trials = 1 << 12;

class Matchfinder
  {
  enum { num_prev_positions = 1 << (8 * min_match_len) };

  long long partial_file_pos;
  const int dictionary_size;	// bytes to keep in buffer before pos
  const int dictionary_mask;
  const int after_size;		// bytes to keep in buffer after pos
  const int buffer_size;
  uint8_t * const buffer;
  int pos;
  int stream_pos;		// first byte not yet read from file
  int pos_limit;		// when reached, a new block must be read
  const int _ides;
  uint32_t _crc;
  int cyclic_pos;
  const int match_len_limit;
  int * const prev_positions;	// last seen position of key
  int * const prev_pos_tree;
  bool at_stream_end;		// stream_pos shows real end of file

  bool read_block() throw();

public:
  Matchfinder( const int dictionary_bits, const int len_limit, const int ides )
    :
    partial_file_pos( 0 ),
    dictionary_size( 1 << dictionary_bits ),
    dictionary_mask( dictionary_size - 1 ),
    after_size( max_match_len ),
    buffer_size( ( 2 * dictionary_size ) + max_num_trials + after_size ),
    buffer( new uint8_t[buffer_size] ),
    pos( 0 ),
    stream_pos( 0 ),
    pos_limit( 0 ),
    _ides( ides ),
    _crc( 0xFFFFFFFF ),
    cyclic_pos( 0 ),
    match_len_limit( len_limit ),
    prev_positions( new int[num_prev_positions] ),
    prev_pos_tree( new int[dictionary_size] ),
    at_stream_end( false )
    {
    if( !read_block() ) throw Error( "read error" );
    for( int i = 0; i < num_prev_positions; ++i ) prev_positions[i] = -1;
    }

  ~Matchfinder() { delete[] prev_pos_tree; delete[] prev_positions; delete[] buffer; }

  uint8_t operator[]( const int i ) const throw() { return buffer[pos+i]; }
  int available_bytes() const throw() { return stream_pos - pos; }
  uint32_t crc() const throw() { return _crc ^ 0xFFFFFFFF; }
  long long file_position() const throw() { return partial_file_pos + pos; }
  bool finished() const throw() { return pos >= stream_pos; }
  bool move_pos() throw();
  const uint8_t * ptr_to_current_pos() const throw() { return buffer + pos; }

  int longest_match_len( int * const distances = 0 ) throw();

  int true_match_len( const int index, const int distance, int len_limit ) const throw()
    {
    if( index + len_limit > available_bytes() )
      len_limit = available_bytes() - index;
    const uint8_t * const data = buffer + pos + index;
    int i = 0;
    while( i < len_limit && data[i] == data[i-distance] ) ++i;
    return i;
    }
  };


const int price_shift = 6;

class Prob_prices
  {
  int data[bit_model_total >> 2];

public:
  Prob_prices()
    {
    const int num_bits = ( bit_model_total_bits - 2 );
    for( int i = num_bits - 1; i >= 0; --i )
      {
      int start = 1 << ( num_bits - i - 1 );
      int end = 1 << ( num_bits - i);
      for( int j = start; j < end; ++j )
        data[j] = (i << price_shift) +
            ( ((end - j) << price_shift) >> (num_bits - i - 1) );
      }
    }

  int operator[]( const int symbol ) const throw()
    { return data[symbol >> 2]; }
  };

extern const Prob_prices prob_prices;

inline int price0( const Bit_model & bm ) throw()
  { return prob_prices[bm.probability]; }

inline int price1( const Bit_model & bm ) throw()
  { return prob_prices[bit_model_total-bm.probability]; }

inline int price_bit( const Bit_model & bm, const int bit ) throw()
  { if( bit ) return price1( bm ); else return price0( bm ); }

inline int price_symbol( const Bit_model bm[], int symbol, const int num_bits ) throw()
  {
  symbol |= ( 1 << num_bits );
  int price = 0;
  while( symbol > 1 )
    {
    price += price_bit( bm[(symbol>>1)-1], symbol & 1 );
    symbol >>= 1;
    }
  return price;
  }

inline int price_symbol_reversed( const Bit_model bm[], int symbol,
                                  const int num_bits ) throw()
  {
  int price = 0;
  int model = 1;
  for( int i = num_bits; i > 0; --i )
    {
    const int bit = symbol & 1;
    symbol >>= 1;
    price += price_bit( bm[model-1], bit );
    model = ( model << 1 ) | bit;
    }
  return price;
  }

inline int price_matched( const Bit_model bm[], const int symbol,
                          const int match_byte ) throw()
  {
  int price = 0;
  int model = 1;

  for( int i = 7; i >= 0; --i )
    {
    const int match_bit = ( match_byte >> i ) & 1;
    const int bit = ( symbol >> i ) & 1;
    price += price_bit( bm[(match_bit<<8)+model+0xFF], bit );
    model = ( model << 1 ) | bit;
    if( match_bit != bit )
      while( --i >= 0 )
      {
      const int bit = ( symbol >> i ) & 1;
      price += price_bit( bm[model-1], bit );
      model = ( model << 1 ) | bit;
      }
    }
  return price;
  }


class Range_encoder
  {
  uint64_t low;
  long long partial_file_pos;
  const int buffer_size;
  uint8_t * const buffer;
  int pos;
  uint32_t range;
  int ff_count;
  const int _odes;
  uint8_t cache;

  void shift_low()
    {
    const uint32_t carry = low >> 32;
    if( low < 0xFF000000 || carry == 1 )
      {
      put_byte( cache + carry );
      for( ; ff_count > 0; --ff_count ) put_byte( 0xFF + carry );
      cache = low >> 24;
      }
    else ++ff_count;
    low = (uint32_t)low << 8;
    }

public:
  Range_encoder( const int header_size, const int odes )
    :
    low( 0 ),
    partial_file_pos( header_size ),
    buffer_size( 65536 ),
    buffer( new uint8_t[buffer_size] ),
    pos( 0 ),
    range( 0xFFFFFFFF ),
    ff_count( 0 ),
    _odes( odes ),
    cache( 0 ) {}

  ~Range_encoder() { delete[] buffer; }

  void put_byte( const uint8_t b )
    {
    buffer[pos] = b;
    if( ++pos >= buffer_size ) flush();
    }

  void encode( const int symbol, const int num_bits = 1 )
    {
    for( int i = num_bits - 1; i >= 0; --i )
      {
      range >>= 1;
      if( (symbol >> i) & 1 ) low += range;
      if( range <= 0xFFFFFF ) { range <<= 8; shift_low(); }
      }
    }

  void encode_bit( Bit_model & bm, const int bit )
    {
    const uint32_t bound = ( range >> bit_model_total_bits ) * bm.probability;
    if( !bit )
      {
      range = bound;
      bm.probability += (bit_model_total - bm.probability) >> bit_model_move_bits;
      }
    else
      {
      low += bound;
      range -= bound;
      bm.probability -= bm.probability >> bit_model_move_bits;
      }
    if( range <= 0xFFFFFF ) { range <<= 8; shift_low(); }
    }

  void encode_tree( Bit_model bm[], const int symbol, const int num_bits )
    {
    int mask = ( 1 << ( num_bits - 1 ) );
    int model = 1;
    for( int i = num_bits; i > 0; --i, mask >>= 1 )
      {
      const int bit = ( symbol & mask );
      encode_bit( bm[model-1], bit );
      model <<= 1;
      if( bit ) model |= 1;
      }
    }

  void encode_tree_reversed( Bit_model bm[], int symbol, const int num_bits )
    {
    int model = 1;
    for( int i = 0; i < num_bits; ++i )
      {
      const int bit = symbol & 1;
      encode_bit( bm[model-1], bit );
      model <<= 1;
      if( bit ) model |= 1;
      symbol >>= 1;
      }
    }

  void encode_matched( Bit_model bm[], int symbol, int match_byte )
    {
    int model = 1;
    for( int i = 7; i >= 0; --i )
      {
      const int bit = ( symbol >> i ) & 1;
      const int match_bit = ( match_byte >> i ) & 1;
      encode_bit( bm[(match_bit<<8)+model+0xFF], bit );
      model = ( model << 1 ) | bit;
      if( match_bit != bit )
        while( --i >= 0 )
          {
          const int bit = ( symbol >> i ) & 1;
          encode_bit( bm[model-1], bit );
          model = ( model << 1 ) | bit;
          }
      }
    }

  void flush()
    {
    if( pos > 0 )
      {
      const int wr = writeblock( _odes, (char *)buffer, pos );
      if( wr != pos ) throw Error( "write error" );
      partial_file_pos += pos;
      pos = 0;
      }
    }

  void flush_data() { for( int i = 0; i < 5; ++i ) shift_low(); }

  long long file_position() const throw()
    { return partial_file_pos + pos + ff_count; }
  };


class Len_encoder
  {
  Bit_model choice1;
  Bit_model choice2;
  Bit_model bm_low[pos_states][len_low_symbols];
  Bit_model bm_mid[pos_states][len_mid_symbols];
  Bit_model bm_high[len_high_symbols];
  int prices[pos_states][max_len_symbols];
  const int len_symbols;
  int counters[pos_states];

  int calculate_price( const int symbol, const int pos_state ) const throw();

  void update_prices( const int pos_state ) throw()
    {
    for( int len = 0; len < len_symbols; ++len )
      prices[pos_state][len] = calculate_price( len, pos_state );
    counters[pos_state] = len_symbols;
    }

public:
  Len_encoder( const int len_limit )
    : len_symbols( len_limit + 1 - min_match_len )
    {
    for( int i = 0; i < pos_states; ++i ) update_prices( i );
    }

  void encode( Range_encoder & range_encoder, const int symbol,
               const int pos_state );

  int price( const int symbol, const int pos_state ) const throw()
    { return prices[pos_state][symbol]; }
  };


class Literal_encoder
  {
  typedef Bit_model Bm_array[0x300];
  Bm_array * const bm_literal;

  int state( const int prev_byte ) const throw()
    { return ( prev_byte >> ( 8 - literal_context_bits ) ); }

public:
  Literal_encoder()
    : bm_literal( new Bm_array[1<<literal_context_bits] ) {}

  ~Literal_encoder() { delete[] bm_literal; }

  void encode( Range_encoder & range_encoder, uint8_t prev_byte, uint8_t symbol )
    { range_encoder.encode_tree( bm_literal[state(prev_byte)], symbol, 8 ); }

  void encode_matched( Range_encoder & range_encoder, uint8_t prev_byte, uint8_t match_byte, uint8_t symbol )
    { range_encoder.encode_matched( bm_literal[state(prev_byte)], symbol, match_byte ); }

  int price_matched( uint8_t prev_byte, uint8_t symbol, uint8_t match_byte ) const throw()
    { return ::price_matched( bm_literal[state(prev_byte)], symbol, match_byte ); }

  int price_symbol( uint8_t prev_byte, uint8_t symbol ) const throw()
    { return ::price_symbol( bm_literal[state(prev_byte)], symbol, 8 ); }
  };


class LZ_encoder
  {
  enum { dis_align_mask = dis_align_size - 1,
         infinite_price = 0x0FFFFFFF,
         num_rep_distances = 4 };	// must be 4

  struct Trial
    {
    State state;
    int dis;
    int prev_index;
    int price;		// dual use var; cumulative price, match length
    int reps[num_rep_distances];
    void update( const int d, const int p_i, const int pr )
      { if( pr < price ) { dis = d; prev_index = p_i; price = pr; } }
    };

  const int match_len_limit;
  const int num_dis_slots;
  int longest_match_found;

  Bit_model bm_match[State::states][pos_states];
  Bit_model bm_rep[State::states];
  Bit_model bm_rep0[State::states];
  Bit_model bm_rep1[State::states];
  Bit_model bm_rep2[State::states];
  Bit_model bm_len[State::states][pos_states];
  Bit_model bm_dis_slot[max_dis_states][1<<dis_slot_bits];
  Bit_model bm_dis[modeled_distances-end_dis_model];
  Bit_model bm_align[dis_align_size];

  Matchfinder matchfinder;
  Range_encoder range_encoder;
  Len_encoder len_encoder;
  Len_encoder rep_match_len_encoder;
  Literal_encoder literal_encoder;

  int dis_slots[1024];
  int match_distances[max_match_len+1];
  Trial trials[max_num_trials];

  int dis_slot_prices[max_dis_states][2*max_dictionary_bits];
  int dis_prices[max_dis_states][modeled_distances];
  int align_prices[dis_align_size];
  int align_price_count;

  void fill_align_prices() throw();
  void fill_distance_prices() throw();

  int get_dis_slot( const unsigned int dis ) const throw()
    {
    if( dis < (1 << 10) ) return dis_slots[dis];
    if( dis < (1 << 19) ) return dis_slots[dis>>9] + 18;
    if( dis < (1 << 28) ) return dis_slots[dis>>18] + 36;
    return dis_slots[dis>>27] + 54;
    }

  void mtf_reps( const int dis, int reps[num_rep_distances] ) throw()
    {
    if( dis >= num_rep_distances )
      {
      for( int i = num_rep_distances - 1; i > 0; --i ) reps[i] = reps[i-1];
      reps[0] = dis - num_rep_distances;
      }
    else if( dis > 0 )
      {
      const int distance = reps[dis];
      for( int i = dis; i > 0; --i ) reps[i] = reps[i-1];
      reps[0] = distance;
      }
    }

  int price_rep_len1( const State & state, const int pos_state ) const throw()
    {
    return price0( bm_rep0[state()] ) + price0( bm_len[state()][pos_state] );
    }

  int price_rep( const int rep, const int len, const State & state,
                 const int pos_state ) const throw()
    {
    int price = rep_match_len_encoder.price( len - min_match_len, pos_state );
    if( rep == 0 )
      {
      price += price0( bm_rep0[state()] );
      price += price1( bm_len[state()][pos_state] );
      }
    else
      {
      price += price1( bm_rep0[state()] );
      if( rep == 1 )
        price += price0( bm_rep1[state()] );
      else
        {
        price += price1( bm_rep1[state()] );
        price += price_bit( bm_rep2[state()], rep - 2 );
        }
      }
    return price;
    }

  int price_pair( const int dis, const int len, const int pos_state ) const throw()
    {
    if( len <= min_match_len && dis >= 256 ) return infinite_price;
    int price;
    const int dis_state = get_dis_state( len );
    if( dis < modeled_distances )
      price = dis_prices[dis_state][dis];
    else
      price = dis_slot_prices[dis_state][get_dis_slot(dis)] +
              align_prices[dis & dis_align_mask];
    return price + len_encoder.price( len - min_match_len, pos_state );
    }

  int read_match_distances() throw()
    {
    int len = matchfinder.longest_match_len( match_distances );
    if( len == match_len_limit )
      len += matchfinder.true_match_len( len, match_distances[len] + 1, max_match_len - len );
    return len;
    }

  bool move_pos( int n ) throw()
    {
    while( --n >= 0 )
      {
      matchfinder.longest_match_len();
      if( !matchfinder.move_pos() ) return false;
      }
    return true;
    }

  void backward( int cur )
    {
    int & dis = trials[cur].dis;
    while( cur > 0 )
      {
      const int prev_index = trials[cur].prev_index;
      Trial & prev_trial = trials[prev_index];
      std::swap( dis, prev_trial.dis );
      prev_trial.price = cur - prev_index;		// len
      cur = prev_index;
      }
    }

  int best_pair_sequence( const int reps[num_rep_distances],
                          const State & state );

  void flush( const State & state );

public:
  LZ_encoder( const File_header & header, const int ides,
              const int odes, const int len_limit );

  bool encode();

  long long input_file_position() const throw()
    { return matchfinder.file_position(); }
  long long output_file_position() const throw()
    { return range_encoder.file_position(); }
  };
