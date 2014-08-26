/*  Lzip - LZMA lossless data compressor
    Copyright (C) 2008-2014 Antonio Diaz Diaz.

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

enum { max_num_trials = 1 << 13,
       price_shift_bits = 6,
       price_step_bits = 2,
       price_step = 1 << price_step_bits };

class Dis_slots
  {
  uint8_t data[1<<10];

public:
  void init()
    {
    for( int slot = 0; slot < 4; ++slot ) data[slot] = slot;
    for( int i = 4, size = 2, slot = 4; slot < 20; slot += 2 )
      {
      std::memset( &data[i], slot, size );
      std::memset( &data[i+size], slot + 1, size );
      size <<= 1;
      i += size;
      }
    }

  uint8_t operator[]( const int dis ) const { return data[dis]; }
  };

extern Dis_slots dis_slots;

inline uint8_t get_slot( const unsigned dis )
  {
  if( dis < (1 << 10) ) return dis_slots[dis];
  if( dis < (1 << 19) ) return dis_slots[dis>> 9] + 18;
  if( dis < (1 << 28) ) return dis_slots[dis>>18] + 36;
  return dis_slots[dis>>27] + 54;
  }


class Prob_prices
  {
  short data[bit_model_total >> price_step_bits];

public:
  void init()
    {
    for( int i = 0; i < bit_model_total >> price_step_bits; ++i )
      {
      unsigned val = ( i * price_step ) + ( price_step / 2 );
      int bits = 0;				// base 2 logarithm of val
      for( int j = 0; j < price_shift_bits; ++j )
        {
        val = val * val;
        bits <<= 1;
        while( val >= 1 << 16 ) { val >>= 1; ++bits; }
        }
      bits += 15;				// remaining bits in val
      data[i] = ( bit_model_total_bits << price_shift_bits ) - bits;
      }
    }

  int operator[]( const int probability ) const
    { return data[probability >> price_step_bits]; }
  };

extern Prob_prices prob_prices;


inline int price0( const Bit_model bm )
  { return prob_prices[bm.probability]; }

inline int price1( const Bit_model bm )
  { return prob_prices[bit_model_total - bm.probability]; }

inline int price_bit( const Bit_model bm, const int bit )
  { if( bit ) return price1( bm ); else return price0( bm ); }


inline int price_symbol3( const Bit_model bm[], int symbol )
  {
  int bit = symbol & 1;
  symbol |= 8; symbol >>= 1;
  int price = price_bit( bm[symbol], bit );
  bit = symbol & 1; symbol >>= 1; price += price_bit( bm[symbol], bit );
  return price + price_bit( bm[1], symbol & 1 );
  }


inline int price_symbol6( const Bit_model bm[], int symbol )
  {
  int bit = symbol & 1;
  symbol |= 64; symbol >>= 1;
  int price = price_bit( bm[symbol], bit );
  bit = symbol & 1; symbol >>= 1; price += price_bit( bm[symbol], bit );
  bit = symbol & 1; symbol >>= 1; price += price_bit( bm[symbol], bit );
  bit = symbol & 1; symbol >>= 1; price += price_bit( bm[symbol], bit );
  bit = symbol & 1; symbol >>= 1; price += price_bit( bm[symbol], bit );
  return price + price_bit( bm[1], symbol & 1 );
  }


inline int price_symbol8( const Bit_model bm[], int symbol )
  {
  int bit = symbol & 1;
  symbol |= 0x100; symbol >>= 1;
  int price = price_bit( bm[symbol], bit );
  bit = symbol & 1; symbol >>= 1; price += price_bit( bm[symbol], bit );
  bit = symbol & 1; symbol >>= 1; price += price_bit( bm[symbol], bit );
  bit = symbol & 1; symbol >>= 1; price += price_bit( bm[symbol], bit );
  bit = symbol & 1; symbol >>= 1; price += price_bit( bm[symbol], bit );
  bit = symbol & 1; symbol >>= 1; price += price_bit( bm[symbol], bit );
  bit = symbol & 1; symbol >>= 1; price += price_bit( bm[symbol], bit );
  return price + price_bit( bm[1], symbol & 1 );
  }


inline int price_symbol_reversed( const Bit_model bm[], int symbol,
                                  const int num_bits )
  {
  int price = 0;
  int model = 1;
  for( int i = num_bits; i > 0; --i )
    {
    const int bit = symbol & 1;
    price += price_bit( bm[model], bit );
    model = ( model << 1 ) | bit;
    symbol >>= 1;
    }
  return price;
  }


inline int price_matched( const Bit_model bm[], int symbol, int match_byte )
  {
  int price = 0;
  int mask = 0x100;
  symbol |= mask;

  do {
    match_byte <<= 1;
    const int match_bit = match_byte & mask;
    symbol <<= 1;
    const int bit = symbol & 0x100;
    price += price_bit( bm[match_bit+(symbol>>9)+mask], bit );
    mask &= ~(match_byte ^ symbol);	// if( match_bit != bit ) mask = 0;
    }
  while( symbol < 0x10000 );
  return price;
  }


class Matchfinder_base
  {
  bool read_block();
  void normalize_pos();

  Matchfinder_base( const Matchfinder_base & );	// declared as private
  void operator=( const Matchfinder_base & );	// declared as private

protected:
  unsigned long long partial_data_pos;
  uint8_t * buffer;		// input buffer
  int32_t * prev_positions;	// 1 + last seen position of key. else 0
  int32_t * pos_array;		// may be tree or chain
  const int before_size;	// bytes to keep in buffer before dictionary
  int buffer_size;
  int dictionary_size_;		// bytes to keep in buffer before pos
  int pos;			// current pos in buffer
  int cyclic_pos;		// cycles through [0, dictionary_size]
  int stream_pos;		// first byte not yet read from file
  int pos_limit;		// when reached, a new block must be read
  int key4_mask;
  int num_prev_positions;	// size of prev_positions
  int pos_array_size;
  const int infd;		// input file descriptor
  bool at_stream_end;		// stream_pos shows real end of file

  Matchfinder_base( const int before, const int dict_size,
                    const int after_size, const int dict_factor,
                    const int num_prev_positions23,
                    const int pos_array_factor, const int ifd );

  ~Matchfinder_base()
    { delete[] prev_positions; std::free( buffer ); }

public:
  uint8_t operator[]( const int distance ) const
    { return buffer[pos-distance]; }
  int available_bytes() const { return stream_pos - pos; }
  unsigned long long data_position() const { return partial_data_pos + pos; }
  int dictionary_size() const { return dictionary_size_; }
  bool finished() const { return at_stream_end && pos >= stream_pos; }
  const uint8_t * ptr_to_current_pos() const { return buffer + pos; }

  int true_match_len( const int index, const int distance, int len_limit ) const
    {
    if( index + len_limit > available_bytes() )
      len_limit = available_bytes() - index;
    const uint8_t * const data = buffer + pos + index;
    int i = 0;
    while( i < len_limit && data[i-distance] == data[i] ) ++i;
    return i;
    }

  void move_pos()
    {
    if( ++cyclic_pos > dictionary_size_ ) cyclic_pos = 0;
    if( ++pos >= pos_limit ) normalize_pos();
    }

  void reset();
  };


struct Pair			/* distance-length pair */
  {
  int dis;
  int len;
  };


class Matchfinder : public Matchfinder_base
  {
  enum { before_size = max_num_trials + 1,
         // bytes to keep in buffer after pos
         after_size = ( 2 * max_match_len ) + 1,
         dict_factor = 2,
         num_prev_positions3 = 1 << 16,
         num_prev_positions2 = 1 << 10,
         num_prev_positions23 = num_prev_positions2 + num_prev_positions3,
         pos_array_factor = 2 };

  const int cycles;
  const int match_len_limit_;

public:
  Matchfinder( const int dict_size, const int len_limit, const int ifd )
    :
    Matchfinder_base( before_size, dict_size, after_size, dict_factor,
                      num_prev_positions23, pos_array_factor, ifd ),
    cycles( ( len_limit < max_match_len ) ? 16 + ( len_limit / 2 ) : 256 ),
    match_len_limit_( len_limit )
    {}

  bool dec_pos( const int ahead )
    {
    if( ahead < 0 || pos < ahead ) return false;
    pos -= ahead;
    cyclic_pos -= ahead;
    if( cyclic_pos < 0 ) cyclic_pos += dictionary_size_ + 1;
    return true;
    }

  int match_len_limit() const { return match_len_limit_; }
  int get_match_pairs( Pair * pairs = 0 );
  };


class Range_encoder
  {
  enum { buffer_size = 65536 };
  uint64_t low;
  unsigned long long partial_member_pos;
  uint8_t * const buffer;	// output buffer
  int pos;			// current pos in buffer
  uint32_t range;
  unsigned ff_count;
  const int outfd;		// output file descriptor
  uint8_t cache;

  void shift_low()
    {
    const bool carry = ( low > 0xFFFFFFFFU );
    if( carry || low < 0xFF000000U )
      {
      put_byte( cache + carry );
      for( ; ff_count > 0; --ff_count ) put_byte( 0xFF + carry );
      cache = low >> 24;
      }
    else ++ff_count;
    low = ( low & 0x00FFFFFFU ) << 8;
    }

  Range_encoder( const Range_encoder & );	// declared as private
  void operator=( const Range_encoder & );	// declared as private

public:
  explicit Range_encoder( const int ofd )
    :
    low( 0 ),
    partial_member_pos( 0 ),
    buffer( new uint8_t[buffer_size] ),
    pos( 0 ),
    range( 0xFFFFFFFFU ),
    ff_count( 0 ),
    outfd( ofd ),
    cache( 0 )
    {}

  ~Range_encoder() { delete[] buffer; }

  unsigned long long member_position() const
    { return partial_member_pos + pos + ff_count; }

  void flush() { for( int i = 0; i < 5; ++i ) shift_low(); }
  void flush_data();

  void put_byte( const uint8_t b )
    {
    buffer[pos] = b;
    if( ++pos >= buffer_size ) flush_data();
    }

  void encode( const int symbol, const int num_bits )
    {
    for( int i = num_bits - 1; i >= 0; --i )
      {
      range >>= 1;
      if( (symbol >> i) & 1 ) low += range;
      if( range <= 0x00FFFFFFU ) { range <<= 8; shift_low(); }
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
    if( range <= 0x00FFFFFFU ) { range <<= 8; shift_low(); }
    }

  void encode_tree3( Bit_model bm[], const int symbol )
    {
    int model = 1;
    int bit = ( symbol >> 2 ) & 1;
    encode_bit( bm[model], bit ); model = ( model << 1 ) | bit;
    bit = ( symbol >> 1 ) & 1;
    encode_bit( bm[model], bit ); model = ( model << 1 ) | bit;
    encode_bit( bm[model], symbol & 1 );
    }

  void encode_tree6( Bit_model bm[], const int symbol )
    {
    int model = 1;
    int bit = ( symbol >> 5 ) & 1;
    encode_bit( bm[model], bit ); model = ( model << 1 ) | bit;
    bit = ( symbol >> 4 ) & 1;
    encode_bit( bm[model], bit ); model = ( model << 1 ) | bit;
    bit = ( symbol >> 3 ) & 1;
    encode_bit( bm[model], bit ); model = ( model << 1 ) | bit;
    bit = ( symbol >> 2 ) & 1;
    encode_bit( bm[model], bit ); model = ( model << 1 ) | bit;
    bit = ( symbol >> 1 ) & 1;
    encode_bit( bm[model], bit ); model = ( model << 1 ) | bit;
    encode_bit( bm[model], symbol & 1 );
    }

  void encode_tree8( Bit_model bm[], const int symbol )
    {
    int model = 1;
    int mask = ( 1 << 7 );
    do {
      const int bit = ( symbol & mask );
      encode_bit( bm[model], bit );
      model <<= 1; if( bit ) ++model;
      }
    while( mask >>= 1 );
    }

  void encode_tree_reversed( Bit_model bm[], int symbol, const int num_bits )
    {
    int model = 1;
    for( int i = num_bits; i > 0; --i )
      {
      const int bit = symbol & 1;
      encode_bit( bm[model], bit );
      model = ( model << 1 ) | bit;
      symbol >>= 1;
      }
    }

  void encode_matched( Bit_model bm[], int symbol, int match_byte )
    {
    int mask = 0x100;
    symbol |= mask;

    do {
      match_byte <<= 1;
      const int match_bit = match_byte & mask;
      symbol <<= 1;
      const int bit = symbol & 0x100;
      encode_bit( bm[match_bit+(symbol>>9)+mask], bit );
      mask &= ~(match_byte ^ symbol);	// if( match_bit != bit ) mask = 0;
      }
    while( symbol < 0x10000 );
    }

  void encode_len( Len_model & lm, int symbol, const int pos_state )
    {
    symbol -= min_match_len;
    bool bit = ( symbol >= len_low_symbols );
    encode_bit( lm.choice1, bit );
    if( !bit )
      encode_tree3( lm.bm_low[pos_state], symbol );
    else
      {
      bit = ( symbol >= len_low_symbols + len_mid_symbols );
      encode_bit( lm.choice2, bit );
      if( !bit )
        encode_tree3( lm.bm_mid[pos_state], symbol - len_low_symbols );
      else
        encode_tree8( lm.bm_high, symbol - len_low_symbols - len_mid_symbols );
      }
    }
  };


class Len_prices
  {
  const Len_model & lm;
  const int len_symbols;
  const int count;
  int prices[pos_states][max_len_symbols];
  int counters[pos_states];

  void update_low_mid_prices( const int pos_state )
    {
    counters[pos_state] = count;
    int * const pps = prices[pos_state];
    int tmp = price0( lm.choice1 );
    int len = 0;
    for( ; len < len_low_symbols && len < len_symbols; ++len )
      pps[len] = tmp + price_symbol3( lm.bm_low[pos_state], len );
    if( len >= len_symbols ) return;
    tmp = price1( lm.choice1 ) + price0( lm.choice2 );
    for( ; len < len_low_symbols + len_mid_symbols && len < len_symbols; ++len )
      pps[len] = tmp +
                 price_symbol3( lm.bm_mid[pos_state], len - len_low_symbols );
    }

  void update_high_prices()
    {
    const int tmp = price1( lm.choice1 ) + price1( lm.choice2 );
    for( int len = len_low_symbols + len_mid_symbols; len < len_symbols; ++len )
      // using 4 slots per value makes "price" faster
      prices[3][len] = prices[2][len] = prices[1][len] = prices[0][len] = tmp +
        price_symbol8( lm.bm_high, len - len_low_symbols - len_mid_symbols );
    }

public:
  Len_prices( const Len_model & m, const int match_len_limit )
    :
    lm( m ),
    len_symbols( match_len_limit + 1 - min_match_len ),
    count( ( match_len_limit > 12 ) ? 1 : len_symbols )
    {
    for( int i = 0; i < pos_states; ++i ) counters[i] = 0;
    }

  void decrement_counter( const int pos_state ) { --counters[pos_state]; }

  void update_prices()
    {
    bool high_pending = false;
    for( int pos_state = 0; pos_state < pos_states; ++pos_state )
      if( counters[pos_state] <= 0 )
        { update_low_mid_prices( pos_state ); high_pending = true; }
    if( high_pending && len_symbols > len_low_symbols + len_mid_symbols )
      update_high_prices();
    }

  int price( const int symbol, const int pos_state ) const
    { return prices[pos_state][symbol - min_match_len]; }
  };


class LZ_encoder_base
  {
protected:
  enum { max_marker_size = 16,
         num_rep_distances = 4 };	// must be 4

  uint32_t crc_;

  Bit_model bm_literal[1<<literal_context_bits][0x300];
  Bit_model bm_match[State::states][pos_states];
  Bit_model bm_rep[State::states];
  Bit_model bm_rep0[State::states];
  Bit_model bm_rep1[State::states];
  Bit_model bm_rep2[State::states];
  Bit_model bm_len[State::states][pos_states];
  Bit_model bm_dis_slot[len_states][1<<dis_slot_bits];
  Bit_model bm_dis[modeled_distances-end_dis_model];
  Bit_model bm_align[dis_align_size];

  Range_encoder renc;
  Len_model match_len_model;
  Len_model rep_len_model;

  unsigned crc() const { return crc_ ^ 0xFFFFFFFFU; }

  LZ_encoder_base( const File_header & header, const int outfd )
    :
    crc_( 0xFFFFFFFFU ),
    renc( outfd )
    {
    for( int i = 0; i < File_header::size; ++i )
      renc.put_byte( header.data[i] );
    }

  int price_literal( const uint8_t prev_byte, const uint8_t symbol ) const
    { return price_symbol8( bm_literal[get_lit_state(prev_byte)], symbol ); }

  int price_matched( const uint8_t prev_byte, const uint8_t symbol,
                     const uint8_t match_byte ) const
    { return ::price_matched( bm_literal[get_lit_state(prev_byte)], symbol,
                              match_byte ); }

  void encode_literal( const uint8_t prev_byte, const uint8_t symbol )
    { renc.encode_tree8( bm_literal[get_lit_state(prev_byte)], symbol ); }

  void encode_matched( const uint8_t prev_byte, const uint8_t symbol,
                       const uint8_t match_byte )
    { renc.encode_matched( bm_literal[get_lit_state(prev_byte)], symbol,
                           match_byte ); }

  void encode_pair( const unsigned dis, const int len, const int pos_state )
    {
    renc.encode_len( match_len_model, len, pos_state );
    const int dis_slot = get_slot( dis );
    renc.encode_tree6( bm_dis_slot[get_len_state(len)], dis_slot );

    if( dis_slot >= start_dis_model )
      {
      const int direct_bits = ( dis_slot >> 1 ) - 1;
      const unsigned base = ( 2 | ( dis_slot & 1 ) ) << direct_bits;
      const unsigned direct_dis = dis - base;

      if( dis_slot < end_dis_model )
        renc.encode_tree_reversed( bm_dis + base - dis_slot - 1, direct_dis,
                                   direct_bits );
      else
        {
        renc.encode( direct_dis >> dis_align_bits, direct_bits - dis_align_bits );
        renc.encode_tree_reversed( bm_align, direct_dis, dis_align_bits );
        }
      }
    }

  void full_flush( const unsigned long long data_position, const State state );

public:
  unsigned long long member_position() const { return renc.member_position(); }
  };


class LZ_encoder : public LZ_encoder_base
  {
  enum { infinite_price = 0x0FFFFFFF,
         single_step_trial = -2,
         dual_step_trial = -1 };

  struct Trial
    {
    State state;
    int price;		// dual use var; cumulative price, match length
    int dis;		// rep index or match distance. (-1 for literal)
    int prev_index;	// index of prev trial in trials[]
    int prev_index2;	//   -2  trial is single step
			//   -1  literal + rep0
			// >= 0  ( rep or match ) + literal + rep0
    int reps[num_rep_distances];

    void update( const int pr, const int distance, const int p_i )
      {
      if( pr < price )
        { price = pr; dis = distance; prev_index = p_i;
          prev_index2 = single_step_trial; }
      }

    void update2( const int pr, const int p_i )
      {
      if( pr < price )
        { price = pr; dis = 0; prev_index = p_i;
          prev_index2 = dual_step_trial; }
      }

    void update3( const int pr, const int distance, const int p_i,
                  const int p_i2 )
      {
      if( pr < price )
        { price = pr; dis = distance; prev_index = p_i; prev_index2 = p_i2; }
      }
    };

  Matchfinder & matchfinder;
  Len_prices match_len_prices;
  Len_prices rep_len_prices;
  int pending_num_pairs;
  Pair pairs[max_match_len+1];
  Trial trials[max_num_trials];

  int dis_slot_prices[len_states][2*max_dictionary_bits];
  int dis_prices[len_states][modeled_distances];
  int align_prices[dis_align_size];
  const int num_dis_slots;

  void update_distance_prices();

       // move-to-front dis in/into reps if( dis > 0 )
  static void mtf_reps( const int dis, int reps[num_rep_distances] )
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

  int price_shortrep( const State state, const int pos_state ) const
    {
    return price0( bm_rep0[state()] ) + price0( bm_len[state()][pos_state] );
    }

  int price_rep( const int rep, const State state, const int pos_state ) const
    {
    if( rep == 0 ) return price0( bm_rep0[state()] ) +
                          price1( bm_len[state()][pos_state] );
    int price = price1( bm_rep0[state()] );
    if( rep == 1 )
      price += price0( bm_rep1[state()] );
    else
      {
      price += price1( bm_rep1[state()] );
      price += price_bit( bm_rep2[state()], rep - 2 );
      }
    return price;
    }

  int price_rep0_len( const int len, const State state, const int pos_state ) const
    {
    return price_rep( 0, state, pos_state ) +
           rep_len_prices.price( len, pos_state );
    }

  int price_pair( const int dis, const int len, const int pos_state ) const
    {
    const int price = match_len_prices.price( len, pos_state );
    const int len_state = get_len_state( len );
    if( dis < modeled_distances )
      return price + dis_prices[len_state][dis];
    else
      return price + dis_slot_prices[len_state][get_slot( dis )] +
             align_prices[dis & (dis_align_size - 1)];
    }

  int read_match_distances()
    {
    const int num_pairs = matchfinder.get_match_pairs( pairs );
    if( num_pairs > 0 )
      {
      int len = pairs[num_pairs-1].len;
      if( len == matchfinder.match_len_limit() && len < max_match_len )
        {
        len += matchfinder.true_match_len( len, pairs[num_pairs-1].dis + 1,
                                           max_match_len - len );
        pairs[num_pairs-1].len = len;
        }
      }
    return num_pairs;
    }

  void move_pos( int n )
    {
    while( true )
      {
      matchfinder.move_pos();
      if( --n <= 0 ) break;
      matchfinder.get_match_pairs();
      }
    }

  void backward( int cur )
    {
    int & dis = trials[cur].dis;
    while( cur > 0 )
      {
      const int prev_index = trials[cur].prev_index;
      Trial & prev_trial = trials[prev_index];

      if( trials[cur].prev_index2 != single_step_trial )
        {
        prev_trial.dis = -1;
        prev_trial.prev_index = prev_index - 1;
        prev_trial.prev_index2 = single_step_trial;
        if( trials[cur].prev_index2 >= 0 )
          {
          Trial & prev_trial2 = trials[prev_index-1];
          prev_trial2.dis = dis; dis = 0;
          prev_trial2.prev_index = trials[cur].prev_index2;
          prev_trial2.prev_index2 = single_step_trial;
          }
        }
      prev_trial.price = cur - prev_index;			// len
      cur = dis; dis = prev_trial.dis; prev_trial.dis = cur;
      cur = prev_index;
      }
    }

  int sequence_optimizer( const int reps[num_rep_distances],
                          const State state );

public:
  LZ_encoder( Matchfinder & mf, const File_header & header, const int outfd )
    :
    LZ_encoder_base( header, outfd ),
    matchfinder( mf ),
    match_len_prices( match_len_model, mf.match_len_limit() ),
    rep_len_prices( rep_len_model, mf.match_len_limit() ),
    pending_num_pairs( 0 ),
    num_dis_slots( 2 * real_bits( mf.dictionary_size() - 1 ) )
    {}

  bool encode_member( const unsigned long long member_size );
  };
