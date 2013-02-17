/*  Lzip - Data compressor based on the LZMA algorithm
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

enum { max_num_trials = 1 << 12,
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

static inline uint8_t get_slot( const uint32_t dis )
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
    for( int i = price_step / 2; i < bit_model_total; i += price_step )
      {
      unsigned val = i;
      int bits = 0;			// base 2 logarithm of val
      for( int j = 0; j < price_shift_bits; ++j )
        {
        val = val * val;
        bits <<= 1;
        while( val >= 1 << 16 ) { val >>= 1; ++bits; }
        }
      bits += 15;			// remaining bits in val
      data[i >> price_step_bits] =
        ( bit_model_total_bits << price_shift_bits ) - bits;
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


inline int price_symbol( const Bit_model bm[], int symbol, const int num_bits )
  {
  int price = 0;
  symbol |= ( 1 << num_bits );
  while( symbol > 1 )
    {
    const int bit = symbol & 1;
    symbol >>= 1;
    price += price_bit( bm[symbol], bit );
    }
  return price;
  }


inline int price_symbol_reversed( const Bit_model bm[], int symbol,
                                  const int num_bits )
  {
  int price = 0;
  int model = 1;
  for( int i = num_bits; i > 0; --i )
    {
    const int bit = symbol & 1;
    symbol >>= 1;
    price += price_bit( bm[model], bit );
    model = ( model << 1 ) | bit;
    }
  return price;
  }


inline int price_matched( const Bit_model bm[], unsigned symbol,
                          unsigned match_byte )
  {
  int price = 0;
  unsigned mask = 0x100;
  symbol |= 0x100;

  do {
    match_byte <<= 1;
    const unsigned match_bit = match_byte & mask;
    symbol <<= 1;
    const unsigned bit = symbol & 0x100;
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
  int32_t * prev_positions;	// last seen position of key
  int32_t * pos_array;		// may be tree or chain
  const int before_size;	// bytes to keep in buffer before dictionary
  const int match_len_limit_;
  int buffer_size;
  int dictionary_size_;		// bytes to keep in buffer before pos
  int pos;			// current pos in buffer
  int cyclic_pos;		// cycles through [0, dictionary_size]
  int pos_limit;		// when reached, a new block must be read
  int stream_pos;		// first byte not yet read from file
  unsigned key4_mask;
  int num_prev_positions;	// size of prev_positions
  int pos_array_size;
  const int infd;		// input file descriptor
  bool at_stream_end;		// stream_pos shows real end of file

  Matchfinder_base( const int before, const int dict_size,
                    const int after_size, const int dict_factor,
                    const int match_len_limit, const int num_prev_positions23,
                    const int pos_array_factor, const int ifd );

  ~Matchfinder_base()
    { delete[] prev_positions; std::free( buffer ); }

public:
  uint8_t operator[]( const int i ) const { return buffer[pos+i]; }
  int available_bytes() const { return stream_pos - pos; }
  unsigned long long data_position() const { return partial_data_pos + pos; }
  int dictionary_size() const { return dictionary_size_; }
  bool finished() const { return at_stream_end && pos >= stream_pos; }
  int match_len_limit() const { return match_len_limit_; }
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

public:
  Matchfinder( const int dict_size, const int match_len_limit, const int ifd )
    :
    Matchfinder_base( before_size, dict_size, after_size, dict_factor,
                      match_len_limit, num_prev_positions23, pos_array_factor, ifd ),
    cycles( ( match_len_limit < max_match_len ) ?
            16 + ( match_len_limit / 2 ) : 256 )
    {}

  bool dec_pos( const int ahead )
    {
    if( ahead < 0 || pos < ahead ) return false;
    pos -= ahead;
    cyclic_pos -= ahead;
    if( cyclic_pos < 0 ) cyclic_pos += dictionary_size_ + 1;
    return true;
    }

  int get_match_pairs( struct Pair * pairs = 0 );
  };


class Range_encoder
  {
  enum { buffer_size = 65536 };
  uint64_t low;
  unsigned long long partial_member_pos;
  uint8_t * const buffer;	// output buffer
  int pos;			// current pos in buffer
  uint32_t range;
  int ff_count;
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

  void encode_tree( Bit_model bm[], const int symbol, const int num_bits )
    {
    int mask = ( 1 << ( num_bits - 1 ) );
    int model = 1;
    for( int i = num_bits; i > 0; --i, mask >>= 1 )
      {
      const int bit = ( symbol & mask );
      encode_bit( bm[model], bit );
      model <<= 1;
      if( bit ) model |= 1;
      }
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

  void encode_matched( Bit_model bm[], unsigned symbol, unsigned match_byte )
    {
    unsigned mask = 0x100;
    symbol |= 0x100;

    do {
      match_byte <<= 1;
      const unsigned match_bit = match_byte & mask;
      symbol <<= 1;
      const unsigned bit = symbol & 0x100;
      encode_bit( bm[match_bit+(symbol>>9)+mask], bit );
      mask &= ~(match_byte ^ symbol);	// if( match_bit != bit ) mask = 0;
      }
    while( symbol < 0x10000 );
    }
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

  void update_prices( const int pos_state )
    {
    int * const pps = prices[pos_state];
    int tmp = price0( choice1 );
    int len = 0;
    for( ; len < len_low_symbols && len < len_symbols; ++len )
      pps[len] = tmp +
                 price_symbol( bm_low[pos_state], len, len_low_bits );
    tmp = price1( choice1 );
    for( ; len < len_low_symbols + len_mid_symbols && len < len_symbols; ++len )
      pps[len] = tmp + price0( choice2 ) +
                 price_symbol( bm_mid[pos_state], len - len_low_symbols, len_mid_bits );
    for( ; len < len_symbols; ++len )
      // using 4 slots per value makes "price" faster
      prices[3][len] = prices[2][len] = prices[1][len] = prices[0][len] =
                 tmp + price1( choice2 ) +
                 price_symbol( bm_high, len - len_low_symbols - len_mid_symbols, len_high_bits );
    counters[pos_state] = len_symbols;
    }

public:
  explicit Len_encoder( const int match_len_limit )
    : len_symbols( match_len_limit + 1 - min_match_len )
    {
    for( int i = 0; i < pos_states; ++i ) update_prices( i );
    }

  void encode( Range_encoder & range_encoder, int symbol,
               const int pos_state );

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
  Bit_model bm_dis_slot[max_dis_states][1<<dis_slot_bits];
  Bit_model bm_dis[modeled_distances-end_dis_model];
  Bit_model bm_align[dis_align_size];

  Range_encoder range_encoder;
  Len_encoder len_encoder;
  Len_encoder rep_match_len_encoder;

  const int num_dis_slots;

  unsigned crc() const { return crc_ ^ 0xFFFFFFFFU; }

  LZ_encoder_base( const File_header & header, const int dictionary_size,
                   const int match_len_limit, const int outfd )
    :
    crc_( 0xFFFFFFFFU ),
    range_encoder( outfd ),
    len_encoder( match_len_limit ),
    rep_match_len_encoder( match_len_limit ),
    num_dis_slots( 2 * real_bits( dictionary_size - 1 ) )
    {
    for( int i = 0; i < File_header::size; ++i )
      range_encoder.put_byte( header.data[i] );
    }

       // move-to-front dis in/into reps
  void mtf_reps( const int dis, int reps[num_rep_distances] )
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

  int price_literal( const uint8_t prev_byte, const uint8_t symbol ) const
    { return price_symbol( bm_literal[get_lit_state(prev_byte)], symbol, 8 ); }

  int price_matched( const uint8_t prev_byte, const uint8_t symbol,
                     const uint8_t match_byte ) const
    { return ::price_matched( bm_literal[get_lit_state(prev_byte)],
                              symbol, match_byte ); }

  void encode_literal( const uint8_t prev_byte, const uint8_t symbol )
    { range_encoder.encode_tree( bm_literal[get_lit_state(prev_byte)], symbol, 8 ); }

  void encode_matched( const uint8_t prev_byte, const uint8_t symbol,
                       const uint8_t match_byte )
    { range_encoder.encode_matched( bm_literal[get_lit_state(prev_byte)],
                                    symbol, match_byte ); }

  void encode_pair( const uint32_t dis, const int len, const int pos_state )
    {
    len_encoder.encode( range_encoder, len, pos_state );
    const int dis_slot = get_slot( dis );
    range_encoder.encode_tree( bm_dis_slot[get_dis_state(len)], dis_slot, dis_slot_bits );

    if( dis_slot >= start_dis_model )
      {
      const int direct_bits = ( dis_slot >> 1 ) - 1;
      const uint32_t base = ( 2 | ( dis_slot & 1 ) ) << direct_bits;
      const uint32_t direct_dis = dis - base;

      if( dis_slot < end_dis_model )
        range_encoder.encode_tree_reversed( bm_dis + base - dis_slot - 1,
                                            direct_dis, direct_bits );
      else
        {
        range_encoder.encode( direct_dis >> dis_align_bits, direct_bits - dis_align_bits );
        range_encoder.encode_tree_reversed( bm_align, direct_dis, dis_align_bits );
        }
      }
    }

  void full_flush( const unsigned long long data_position, const State state );

public:
  unsigned long long member_position() const
    { return range_encoder.member_position(); }
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
    int dis;		// rep index or match distance
    int prev_index;	// index of prev trial in trials[]
    int dis2;
    int prev_index2;	//   -2  trial is single step
			//   -1  literal + rep0
			// >= 0  rep or match + literal + rep0
    int reps[num_rep_distances];

    void update( const int pr, const int d, const int p_i )
      {
      if( pr < price )
        { price = pr; dis = d; prev_index = p_i;
          prev_index2 = single_step_trial; }
      }

    void update2( const int pr, const int d, const int p_i )
      {
      if( pr < price )
        { price = pr; dis = d; prev_index = p_i;
          prev_index2 = dual_step_trial; }
      }

    void update3( const int pr, const int d, const int p_i,
                                const int d2, const int p_i2 )
      {
      if( pr < price )
        { price = pr; dis = d; prev_index = p_i;
          dis2 = d2; prev_index2 = p_i2; }
      }
    };

  Matchfinder & matchfinder;
  int pending_num_pairs;
  struct Pair pairs[max_match_len+1];
  Trial trials[max_num_trials];

  int dis_slot_prices[max_dis_states][2*max_dictionary_bits];
  int dis_prices[max_dis_states][modeled_distances];
  int align_prices[dis_align_size];
  int align_price_count;

  void fill_align_prices();
  void fill_distance_prices();

  int price_rep_len1( const State state, const int pos_state ) const
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
           rep_match_len_encoder.price( len, pos_state );
    }

  int price_dis( const int dis, const int dis_state ) const
    {
    if( dis < modeled_distances )
      return dis_prices[dis_state][dis];
    else
      return dis_slot_prices[dis_state][get_slot( dis )] +
             align_prices[dis & (dis_align_size - 1)];
    }

  int price_pair( const int dis, const int len, const int pos_state ) const
    {
    return len_encoder.price( len, pos_state ) +
           price_dis( dis, get_dis_state( len ) );
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
    if( --n >= 0 ) matchfinder.move_pos();
    while( --n >= 0 )
      {
      matchfinder.get_match_pairs();
      matchfinder.move_pos();
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
          prev_trial2.dis = trials[cur].dis2;
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
    LZ_encoder_base( header, mf.dictionary_size(), mf.match_len_limit(), outfd ),
    matchfinder( mf ),
    pending_num_pairs( 0 ),
    align_price_count( 0 )
    {}

  bool encode_member( const unsigned long long member_size );
  };
